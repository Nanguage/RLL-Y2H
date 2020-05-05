extern crate bio;
extern crate clap;

use std::fs::File;
use std::str;
use std::collections::HashMap;
use std::mem;

use clap::{Arg, App};
use bio::alignment::pairwise::*;
use bio::io::fastq;

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


fn main() {
    let matches = App::new("PET Extract")
        .arg(Arg::with_name("fq1")
             .required(true)
             .help("Fastq file of reads 1."))
        .arg(Arg::with_name("linker")
             .short("l")
             .long("linker")
             .required(true)
             .takes_value(true)
             .help("The linker sequence."))
        .arg(Arg::with_name("flanking")
             .short("f")
             .long("flanking")
             .help("Flanking length."))
        .get_matches();

    let fq1_path = matches.value_of("fq1").unwrap();
    let fq1_file = File::open(fq1_path).unwrap();
    let fq1 = fastq::Reader::new(fq1_file);
    let linker = matches.value_of("linker").unwrap();
    let flanking = matches.value_of("flanking").unwrap_or("13");
    let flanking: u8 = flanking.parse().unwrap();

    let mut records = fq1.records();
    let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};

    let mut freq: HashMap<(u64, u64), u64> = HashMap::new();

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

        // align linker to read
        let mut aligner = Aligner::with_capacity(seq.len(), linker.len(), -1, -1, &score);
        let alignment = aligner.semiglobal(&linker.as_bytes(), &seq);

        // filter out non matched reads
        if (alignment.score as f32) < linker.len() as f32 * 0.6 { continue }
        // filter out incomplete flanking
        if (alignment.ystart as u8) < flanking { continue }
        let s = alignment.ystart - flanking as usize;
        let left = &seq[s..alignment.ystart].to_vec();
        let e = alignment.yend + flanking as usize;
        if e > alignment.ylen { continue }
        let right = &seq[alignment.yend..e].to_vec();
        // count left-right pair
        let mut key: (u64, u64) = (compress_seq(left).unwrap(), compress_seq(right).unwrap());
        if key.0 > key.1 { mem::swap(&mut key.0, &mut key.1) };
        *freq.entry(key).or_insert(0) += 1;
    }


    for (k, v) in freq {
        let seq1 = recover_seq(k.0, flanking);
        let seq2 = recover_seq(k.1, flanking);
        println!("{}\t{}\t{}", seq1, seq2, v);
    }
}
