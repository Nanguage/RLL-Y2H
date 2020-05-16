use std::fs::File;
use std::io::Write;
use std::str;
use std::collections::HashMap;
use std::mem;
use std::fmt;
use std::thread;
use std::sync::mpsc;
use std::sync::{Mutex, Arc};
use std::time::Duration;
use std::io::Read;
use std::collections::HashSet;


extern crate bio;
extern crate clap;
extern crate flate2;
extern crate log;
extern crate simple_logger;

use clap::{Arg, App};
use bio::alignment::pairwise::Aligner;
use bio::alignment::Alignment;
use bio::io::fastq;
use bio::alphabets::dna::revcomp;
use log::{info};
use flate2::read::GzDecoder;


const SEQ_NT4_TABLE: [u64; 256] = [
    0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
];

const IDX_TABLE: [u8; 4] = [
    b'A', b'C', b'G', b'T'
];

fn compress_seq(seq: &[u8]) -> Result<u64, &str> {
    let mut res: u64 = 0;
    let mut res_rc: u64 = 0;
    let end = seq.len() - 1;
    for i in 0..seq.len() {
        if i >= 32 {
            return Err("Seq can't longer than 32.")
        }
        let m = SEQ_NT4_TABLE[seq[i] as usize];
        res |= m << i*2;
        res_rc |= (3 - m) << (end - i)*2;
    }
    if res > res_rc { mem::swap(&mut res, &mut res_rc) };
    Ok(res)
}

fn recover_seq(code: u64, k: u8) -> String {
    let mut chars: Vec<u8> = Vec::with_capacity(k as usize);
    for i in 0..k {
        let mask: u64 = 3 << (i*2);
        let idx = (code & mask) >> (i*2);
        let b = IDX_TABLE[idx as usize];
        chars.push(b);
    }
    String::from_utf8(chars).unwrap()
}


enum ExtractRes {
    Ok(String, String),
    ScoreTooLow,
    LeftTooShort,
    RightTooShort,
}


fn extract_pet(seq: &[u8], pattern: &[u8], flanking: u8, score_ratio_thresh: f32) -> (ExtractRes, Alignment) {
    // align linker to read
    let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
    let mut aligner = Aligner::with_capacity(seq.len(), pattern.len(), -1, -1, score);
    let alignment = aligner.semiglobal(pattern, seq);

    // filter out non matched reads
    if (alignment.score as f32) < pattern.len() as f32 * score_ratio_thresh { 
        return (ExtractRes::ScoreTooLow, alignment)
    }
    // filter out incomplete flanking
    if (alignment.ystart as u8) < flanking {
        return (ExtractRes::LeftTooShort, alignment)
    }
    let s = alignment.ystart - flanking as usize;
    let left = String::from_utf8(seq[s..alignment.ystart].to_vec()).unwrap();
    let e = alignment.yend + flanking as usize;
    if e > alignment.ylen {
        return (ExtractRes::RightTooShort, alignment)
    }
    let right = String::from_utf8(seq[alignment.yend..e].to_vec()).unwrap();

    (ExtractRes::Ok(left, right), alignment)
}


struct ResCounter {
    linker_reads: u64,
    score_too_low: u64,
    left_too_short: u64,
    right_too_short: u64,
}

impl ResCounter {
    fn new() -> Self {
        Self {
            linker_reads: 0,
            score_too_low: 0,
            left_too_short: 0,
            right_too_short: 0,
        }
    }

    fn count(&mut self, res: &ExtractRes) {
        match res{
            ExtractRes::Ok(_, _) =>{ self.linker_reads += 1 },
            ExtractRes::ScoreTooLow =>{ self.score_too_low += 1 },
            ExtractRes::LeftTooShort =>{ self.left_too_short += 1 },
            ExtractRes::RightTooShort =>{ self.right_too_short += 1 },
        }
    }
}

impl fmt::Display for ResCounter {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let total = self.linker_reads + self.score_too_low +
                    self.left_too_short + self.right_too_short;
        let ratio = |c| {
            if total == 0 { return format!("0%"); }
            format!("{:.2}%", ((c*100) as f64) / (total as f64))
        };
        write!(f,
            "Count result:
    linker reads\t{}\t{}
    score too low\t{}\t{}
    left too short\t{}\t{}
    right too short\t{}\t{}
total reads: {}\n",
            self.linker_reads, ratio(self.linker_reads),
            self.score_too_low, ratio(self.score_too_low),
            self.left_too_short, ratio(self.left_too_short),
            self.right_too_short, ratio(self.right_too_short),
            total,
        )
    }
}


fn main() {
    simple_logger::init().unwrap();

    let matches = App::new("PET Extract")
        .arg(Arg::with_name("fq")
             .required(true)
             .help("Fastq file of reads 1."))
        .arg(Arg::with_name("linker")
             .short("l")
             .long("linker")
             .required(true)
             .takes_value(true)
             .help("The linker sequence(Not incluede enzyme)."))
        .arg(Arg::with_name("enzyme")
             .short("e")
             .long("enzyme")
             .required(true)
             .takes_value(true)
             .help("Enzyme recognize site.")
            )
        .arg(Arg::with_name("output_prefix")
             .short("o")
             .long("output_prefix")
             .required(true)
             .takes_value(true)
             .help("Prefix of output files."))
        .arg(Arg::with_name("flanking")
             .short("f")
             .long("flanking")
             .takes_value(true)
             .help("Flanking length."))
        .arg(Arg::with_name("score_ratio_thresh")
             .short("s")
             .long("score_ratio_thresh")
             .takes_value(true)
             .help("Threshold of (align score / pattern length)"))
        .arg(Arg::with_name("align_detail")
             .short("d")
             .long("detail")
             .takes_value(true)
             .help("Output the align detail."))
        .arg(Arg::with_name("threads")
             .short("t")
             .long("threads")
             .takes_value(true)
             .help("Number of threads used for processing reads."))
        .arg(Arg::with_name("wait_timeout")
             .long("wait_timeout")
             .takes_value(true)
             .help("Wait time for end channel timeout."))
        .get_matches();

    let fq_path = matches.value_of("fq").unwrap();
    let out_prefix = matches.value_of("output_prefix").unwrap();
    let linker = matches.value_of("linker").unwrap();
    let enzyme = matches.value_of("enzyme").unwrap_or("GTTGGA");
    let flanking = matches.value_of("flanking").unwrap_or("13");
    let flanking: u8 = flanking.parse().unwrap();
    let score_ratio_thresh = matches.value_of("score_ratio_thresh").unwrap_or("0.6");
    let score_ratio_thresh: f32 = score_ratio_thresh.parse().unwrap();
    let threads = matches.value_of("threads").unwrap_or("1");
    let threads: u8 = threads.parse().unwrap();
    let wait_t = matches.value_of("wait_timeout").unwrap_or("500");
    let wait_t: u64 = wait_t.parse().unwrap();

    let mut detail_file = match matches.value_of("align_detail") {
        Some(p) => Some(File::create(p).unwrap()),
        None => None,
    };

    let fq_file: Box<dyn Read + Send + Sync> = if fq_path.ends_with(".gz") {
        Box::new(GzDecoder::new(File::open(fq_path).unwrap()))
    } else {
        Box::new(File::open(fq_path).unwrap())
    };
    let fq = fastq::Reader::new(fq_file);
    let records = fq.records();

    let mut freq: HashMap<(u64, u64), u64> = HashMap::new();

    let l_vec = linker.as_bytes().to_vec();
    let e_vec = enzyme.as_bytes().to_vec();
    let e_rc = revcomp(&e_vec);
    let l_rc = revcomp(&l_vec);
    let patterns = [
        [&e_vec[..], &l_vec[..], &e_rc].concat(),
        [&e_vec[..], &l_rc[..],  &e_rc].concat(),
    ];
    info!("patterns:\n    {}\n    {}",
        str::from_utf8(&patterns[0]).unwrap(),
        str::from_utf8(&patterns[1]).unwrap(),
    );

    let mut counter = ResCounter::new();
    let records = Arc::new(Mutex::new(records));
    let patterns = Arc::new(patterns);
    let mut handles = vec![];
    let (tx, rx) = mpsc::channel();

    info!("Run with {} threads.", threads);

    for _ in 0..threads {
        let records = Arc::clone(&records);
        let patterns = Arc::clone(&patterns);
        let tx1 = mpsc::Sender::clone(&tx);

        let handle = thread::spawn(move || {
            loop {
                // read seq from fq file
                let rec = {
                    let mut records = records.lock().unwrap();
                    match records.next() {
                        Some(r) => match r {
                            Ok(r_) => r_,
                            Err(e) => panic!("{:?}", e),
                        },
                        None => break
                    }
                };
                let seq = String::from_utf8(rec.seq().to_vec()).unwrap();

                let mut align_res: Vec<(ExtractRes, Alignment)> = Vec::with_capacity(2);
                for pattern in patterns.iter() {
                    align_res.push(extract_pet(seq.as_bytes(), &pattern, flanking, score_ratio_thresh));
                    let res = &align_res[align_res.len()-1].0;
                    match res {
                        ExtractRes::Ok(_, _) => {
                            break
                        },
                        _ => {
                            continue    
                        },
                    }
                }
                let rec_id = String::from(rec.id());
                tx1.send((align_res, rec_id)).unwrap();
            }
        });
        handles.push(handle);
    }

    loop {
        match rx.recv_timeout(Duration::from_millis(wait_t)) {
            Ok((align_res, rec_id)) => {
                let res = &align_res[align_res.len()-1];
                let alignment = &res.1;
                if let Some(mut f) = detail_file {
                    // write align detail
                    let _ = writeln!(f,
                        "{}\t{}\t{}\t{}\t{}",
                        rec_id, align_res.len(),
                        alignment.score, alignment.ystart, alignment.yend,
                    );
                    detail_file = Some(f);
                }

                // count left-right pair
                if let ExtractRes::Ok(left, right) = &res.0 {
                    let mut key: (u64, u64) = (compress_seq(left.as_bytes()).unwrap(),
                                               compress_seq(right.as_bytes()).unwrap());
                    if key.0 > key.1 { mem::swap(&mut key.0, &mut key.1) };
                    *freq.entry(key).or_insert(0) += 1;
                }
                counter.count(&res.0);
            },
            _ => {
                info!("End processing.");
                break;
            }
        }
    }

    for handle in handles {  // wait all threads fishish
        handle.join().unwrap();
    }

    let cnt_path = format!("{}.cnt", out_prefix);
    let fq_out_path = format!("{}.cnt.fq", out_prefix);
    let mut cnt_file = File::create(cnt_path.clone()).unwrap();
    let fq_out_file = File::create(fq_out_path.clone()).unwrap();
    let mut fq_out = fastq::Writer::new(fq_out_file);

    let mut key_set = HashSet::new();
    let mut kv_vec = vec![];
    for (k, v) in &freq {
        key_set.insert(k.0);
        key_set.insert(k.1);
        kv_vec.push((k.0, k.1, v));
    }
    info!("Totally {} kinds of pairs and {} kinds of sequences were founded.", freq.len(), key_set.len());
    kv_vec.sort_by(|a, b| b.2.cmp(a.2));
    info!("Write pair counts to tsv file: {}", cnt_path);
    for (k0, k1, v) in kv_vec {
        let _ = writeln!(cnt_file, "{}\t{}\t{}", k0, k1, v);
    }

    info!("Write sequences to fastq file: {}", fq_out_path);
    let mut key_vec = key_set.into_iter().collect::<Vec<u64>>();
    key_vec.sort();
    for k in key_vec {
        let seq = recover_seq(k, flanking);
        let id = format!("{}", k);
        let qual = vec![b'~'; seq.len()];
        let _ = fq_out.write(
            &id,
            Option::None,
            seq.as_bytes(),
            &qual,
        );
    }

    info!("{}", counter);
}
