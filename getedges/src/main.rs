use std::fs::File;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::str::FromStr;
use std::fmt;

extern crate clap;
extern crate log;
extern crate simple_logger;

use clap::{Arg, App};
use log::{info};


enum Node {
    Bait(String),
    Prey(String),
    NotValid(NotValidType),
}

impl fmt::Display for Node {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let node_str = match self {
            Node::Bait(name) => format!("Bait:{}", name),
            Node::Prey(name) => format!("Prey:{}", name),
            Node::NotValid(nvtp) => match nvtp {
                NotValidType::NotFound => format!("NotFound"),
                NotValidType::MapqTooSmall(s) => format!("MAPQTooSmall:{}", s),
                NotValidType::TooManyMisMatch(s) => format!("TooManyMisMatch:{}", s),
                NotValidType::TooManyAligned(s) => format!("TooManyAligned:{}", s),
            },
        };
        write!(f, "{}", node_str)
    }
}

enum NotValidType {
    NotFound,
    MapqTooSmall(u8),
    TooManyMisMatch(u8),
    TooManyAligned(u8),
}

struct SamRec<'a> {
    qname: u64,
    rname: &'a str,
    mapq: u8,
    n_mismatch: u8,
    n_aligned: u8,
}


fn parse_bwa_sam_rec<'a>(line: &'a String) -> SamRec<'a> {
    let items: Vec<&str> = line.split("\t").collect();
    let qname: u64 = items[0].parse().unwrap();
    let rname = items[2];
    let mapq: u8 = items[4].parse().unwrap();
    let mut nm: u8 = 0;
    let mut na: u8 = if rname == "*" {0} else {1};
    if items.len() >= 12 {
        for i in 12..items.len() {
            if items[i].starts_with("NM") {
                let fields: Vec<&str> = items[i].split(":").collect();
                nm = fields[2].parse().unwrap();
                continue
            }
            if items[i].starts_with("XA") {
                let fields: Vec<&str> = items[i].split(";").collect();
                na += fields.len() as u8;
                break
            }
        }
    }
    SamRec {
        qname: qname,
        rname: rname,
        mapq: mapq,
        n_mismatch: nm,
        n_aligned: na,
    }
}


fn load_sam(path: &str, th_mapq: u8, th_mismatch: u8, th_aligned: u8) -> HashMap<u64, Node> {
    let mut key2node = HashMap::new();
    let f = File::open(path).unwrap();
    let buffered = BufReader::new(f);
    for line in buffered.lines() {
        let line = line.unwrap();
        if line.starts_with("@") { continue }
        let rec = parse_bwa_sam_rec(&line);
        let node = if rec.rname == "*" {
            Node::NotValid(NotValidType::NotFound)
        } else if rec.mapq < th_mapq {
            Node::NotValid(NotValidType::MapqTooSmall(rec.mapq))
        } else if rec.n_mismatch > th_mismatch {
            Node::NotValid(NotValidType::TooManyMisMatch(rec.n_mismatch))
        } else if rec.n_aligned > th_aligned {
            Node::NotValid(NotValidType::TooManyAligned(rec.n_aligned))
        } else {
            let name = String::from_str(rec.rname).unwrap();
            if name.starts_with("bait_") {
                Node::Bait(name)
            } else if name.starts_with("prey_") {
                Node::Prey(name)
            } else {
                panic!("Gene name should starts with 'prey_' or 'bait_'.")
            }
        };
        key2node.insert(rec.qname, node);
    }
    key2node
}


struct ResCounter {
    n_valid_pair: u64,
    n_prey_nv_pair: u64,
    n_bait_nv_pair: u64,
    n_nv_pair: u64,
}

impl<'a> ResCounter {
    fn new() -> Self {
        Self {
            n_valid_pair: 0,
            n_prey_nv_pair: 0,
            n_bait_nv_pair: 0,
            n_nv_pair: 0,
        }
    }

    fn count(&mut self,
             bait_prey_cnt: &mut HashMap<(&'a String, &'a String), u64>,
             node1: &'a Node, node2: &'a Node, cnt: u64) {
        match (node1, node2) {
            (Node::Prey(p_name), Node::Bait(b_name)) | (Node::Bait(b_name), Node::Prey(p_name)) => {
                *bait_prey_cnt.entry((b_name, p_name)).or_insert(0) += cnt;
                self.n_valid_pair += 1;
            },
            (Node::Prey(_), Node::NotValid(_)) | (Node::NotValid(_), Node::Prey(_)) => {
                self.n_prey_nv_pair += 1;
            },
            (Node::Bait(_), Node::NotValid(_)) | (Node::NotValid(_), Node::Bait(_)) => {
                self.n_bait_nv_pair += 1;
            },
            (Node::NotValid(_), Node::NotValid(_)) => {
                self.n_nv_pair += 1;
            },
            _ => {}
        }
    }
}

impl fmt::Display for ResCounter {

    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let total = self.n_valid_pair + self.n_bait_nv_pair +
                    self.n_prey_nv_pair + self.n_nv_pair;
        let ratio = |c| {
            if total == 0 { return format!("0%"); }
            format!("{:.2}%", ((c*100) as f64) / (total as f64))
        };
        write!(f,
            "Count result:
    Bait-Prey\t{}\t{}
    Bait-NotValid\t{}\t{}
    Prey-NotValid\t{}\t{}
    NotValid-NotValid\t{}\t{}
total pairs: {}\n",
            self.n_valid_pair, ratio(self.n_valid_pair),
            self.n_bait_nv_pair, ratio(self.n_bait_nv_pair),
            self.n_prey_nv_pair, ratio(self.n_prey_nv_pair),
            self.n_nv_pair, ratio(self.n_nv_pair),
            total,
        )
    }

}


fn main() {
    simple_logger::init().unwrap();
    let matches = App::new("getedges")
        .arg(Arg::with_name("cnt")
             .required(true)
             .help("Pair count file(.cnt)."))
        .arg(Arg::with_name("sam")
             .required(true)
             .help("Aligned sam file"))
        .arg(Arg::with_name("output")
             .short("o")
             .long("output")
             .required(true)
             .takes_value(true)
             .help("Path to output Bait-Prey count pairs, in TSV format."))
        .arg(Arg::with_name("detail")
             .short("d")
             .long("detail")
             .takes_value(true)
             .help("Path to parse detail."))
        .arg(Arg::with_name("th_mapq")
             .long("th_mapq")
             .takes_value(true)
             .help("MAPQ threshold."))
        .arg(Arg::with_name("th_mismatch")
             .long("th_mismatch")
             .takes_value(true)
             .help("Threshold of number of mismatches."))
        .arg(Arg::with_name("th_aligned")
             .long("th_aligned")
             .takes_value(true)
             .help("Threshold of number of alignments."))
        .get_matches();

    let path_cnt = matches.value_of("cnt").unwrap();
    let path_sam = matches.value_of("sam").unwrap();
    let path_out = matches.value_of("output").unwrap();
    let th_mapq: u8 = matches.value_of("th_mapq").unwrap_or("0").parse().unwrap();
    let th_mismatch: u8 = matches.value_of("th_mismatch").unwrap_or("0").parse().unwrap();
    let th_aligned: u8 = matches.value_of("th_aligned").unwrap_or("1").parse().unwrap();
    let mut detail_file = match matches.value_of("detail") {
        Some(p) => Some(File::create(p).unwrap()),
        None => None,
    };

    let key2name = load_sam(path_sam, th_mapq, th_mismatch, th_aligned);
    let cnt_file = BufReader::new(File::open(path_cnt).unwrap());

    let mut bait_prey_cnt:HashMap<(&String, &String), u64> =  HashMap::new();
    let mut res_counter = ResCounter::new();

    for line in cnt_file.lines() {
        let line = line.unwrap();
        let items: Vec<&str> = line.trim().split("\t").collect();
        let key1: u64 = items[0].parse().unwrap();
        let key2: u64 = items[1].parse().unwrap();
        let cnt: u64 = items[2].parse().unwrap();
        let node1 = &key2name[&key1];
        let node2 = &key2name[&key2];

        res_counter.count(&mut bait_prey_cnt, node1, node2, cnt);

        if let Some(mut f) = detail_file {
            // write detail
            let _ = writeln!(f,
                "{}\t{}\t{}\t{}\t{}",
                key1, key2, cnt,
                node1, node2,
            );
            detail_file = Some(f);
        }
    }

    info!("{}", res_counter);

    let mut bait_prey_vec = vec![];
    for ((bait, prey), cnt) in bait_prey_cnt.iter() {
        bait_prey_vec.push(((bait, prey), cnt));
    }
    bait_prey_vec.sort_by(|a, b| b.1.cmp(&a.1));

    let mut file_out = File::create(path_out.clone()).unwrap();
    info!("Output counted Bait-Prey pairs to: {}", path_out);
    for ((bait, prey), cnt) in bait_prey_vec {
        let _ = write!(file_out, "{}\t{}\t{}\n", bait, prey, cnt);
    }

}
