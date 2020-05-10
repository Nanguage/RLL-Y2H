use std::fs::File;
use std::io::Write;
use std::str;
use std::collections::HashMap;
use std::mem;
use std::fmt;

extern crate bio;
extern crate clap;
extern crate simple_logger;

use clap::{Arg, App};
use bio::alignment::pairwise::Aligner;
use bio::alignment::Alignment;
use bio::io::fastq;
use bio::alphabets::dna::revcomp;
use log::{info};


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
    let mut mask: u64;
    for i in 0..seq.len() {
        if i >= 32 {
            return Err("Seq can't longer than 32.")
        }
        mask = SEQ_NT4_TABLE[seq[i] as usize] << i*2;
        res |= mask;
    }
    Ok(res)
}

fn recover_seq(code: u64, k: u8) -> String {
    let mut chars: Vec<u8> = Vec::with_capacity(k as usize);
    for i in 0..(k-1) {
        let mask: u64 = 3 << (i*2);
        let idx = (code & mask) >> (i*2);
        let b = IDX_TABLE[idx as usize];
        chars.push(b);
    }
    String::from_utf8(chars).unwrap()
}


enum ExtractRes<'a> {
    Ok(&'a [u8], &'a [u8]),
    ScoreTooLow,
    LeftTooShort,
    RightTooShort,
}


fn extract_pet<'a>(seq: &'a [u8], pattern: &[u8], flanking: u8) -> (ExtractRes<'a>, Alignment) {
    // align linker to read
    let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
    let mut aligner = Aligner::with_capacity(seq.len(), pattern.len(), -1, -1, score);
    let alignment = aligner.semiglobal(pattern, seq);

    // filter out non matched reads
    if (alignment.score as f32) < pattern.len() as f32 * 0.6 { 
        return (ExtractRes::ScoreTooLow, alignment)
    }
    // filter out incomplete flanking
    if (alignment.ystart as u8) < flanking {
        return (ExtractRes::LeftTooShort, alignment)
    }
    let s = alignment.ystart - flanking as usize;
    let left = &seq[s..alignment.ystart];
    let e = alignment.yend + flanking as usize;
    if e > alignment.ylen {
        return (ExtractRes::RightTooShort, alignment)
    }
    let right = &seq[alignment.yend..e];

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
        .arg(Arg::with_name("output")
             .short("o")
             .long("output")
             .required(true)
             .takes_value(true)
             .help("Output tsv file."))
        .arg(Arg::with_name("flanking")
             .short("f")
             .long("flanking")
             .takes_value(true)
             .help("Flanking length."))
        .arg(Arg::with_name("align_detail")
             .short("d")
             .long("detail")
             .takes_value(true)
             .help("Output the align detail."))
        .get_matches();

    let fq_path = matches.value_of("fq").unwrap();
    let out_path = matches.value_of("output").unwrap();
    let linker = matches.value_of("linker").unwrap();
    let enzyme = matches.value_of("enzyme").unwrap_or("GTTGGA");
    let flanking = matches.value_of("flanking").unwrap_or("13");
    let flanking: u8 = flanking.parse().unwrap();

    let mut detail_file = match matches.value_of("align_detail") {
        Some(p) => Some(File::create(p).unwrap()),
        None => None,
    };

    let fq_file = File::open(fq_path).unwrap();
    let fq = fastq::Reader::new(fq_file);
    let mut out_file = File::create(out_path).unwrap();
    let mut records = fq.records();

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

    loop {
        // read seq from fq file
        let rec = match records.next() {
            Some(r) => match r {
                Ok(r_) => r_,
                Err(e) => panic!("{:?}", e),
            },
            None => break
        };
        let vec = rec.seq();
        let seq = &vec[..(vec.len()-1)];

        let mut align_res: Vec<(ExtractRes, Alignment)> = Vec::with_capacity(2);
        for pattern in patterns.iter() {
            align_res.push(extract_pet(seq, &pattern, flanking));
            let res = &align_res[align_res.len()-1].0;
            match res {
                ExtractRes::Ok(left, right) => {
                    // count left-right pair
                    let mut key: (u64, u64) = (compress_seq(left).unwrap(), compress_seq(right).unwrap());
                    if key.0 > key.1 { mem::swap(&mut key.0, &mut key.1) };
                    *freq.entry(key).or_insert(0) += 1;
                    break
                },
                _ => {
                    continue    
                },
            }
        }
        let alignment = &align_res[align_res.len()-1].1;

        if let Some(mut f) = detail_file {
            // write align detail
            let _ = writeln!(f,
                "{}\t{}\t{}\t{}\t{}",
                rec.id(), align_res.len(),
                alignment.score, alignment.ystart, alignment.yend,
            );
            detail_file = Some(f);
        }

        // count
        counter.count(&align_res[align_res.len()-1].0)
    }

    for (k, v) in freq {
        let seq1 = recover_seq(k.0, flanking);
        let seq2 = recover_seq(k.1, flanking);
        let _ = writeln!(out_file, "{}\t{}\t{}", seq1, seq2, v);
    }

    info!("{}", counter);
}
