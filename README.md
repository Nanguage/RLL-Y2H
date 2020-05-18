# RLL-Y2H

Tools for processing RLL-Y2H data. 

## Install

Download lastest release:

```bash
$ wget http://github.com/Nanguage/RLL-Y2H/releases/download/0.0.1/release.zip
```

Or compile from source code:

```bash
# install rust and cargo
$ curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
# clone repo
$ git clone https://github.com/Nanguage/RLL-Y2H.git
# build from source
$ cd RLL-Y2H
$ cargo build --release
```

## Workflow

### 1. Extract and counting seq pairs

Using the tool 'paircnts'. For example:

```bash
./paircnt ./data/test_R1.fq.gz -l TAGCGTGCGGGTGCCAGGGCGTGCCCTTGAGTTCTCTCAGTTGGGGGCGTTGAC -o test1 -e GTTGGA --threads 10 --flanking 15
```

This will produce two file: `test1.cnt` and `test1.cnt.fq`. The `.cnt` file recorded the 
count of all seq pairs in TSV format. The `.cnt.fq` file is all sequences occured in
the `.cnt` file for find coresponding gene in library by sequence aligment.

More usage detail see:

```bash
./paircnt -h
```


### 2. Make library index and align sequences to the library by bwa aligner

Firstly, make sure bwa aligner is installed on your system. Recormand install with conda:

```bash
$ conda install -c bioconda bwa
```

The library sequencing should stored in `fasta` file, and bait gene's name should starts with
`bait_`, prey gene's name starts with `prey_`, for example:

```
> bait_gene1
aaagcctgcgcatttaattaa
> bait_gene2
attaactgcgcccccaattaa
> prey_gene1
attaactgcgtttttaattga
> prey_gene2
attaagaatccccccaattcc
```

Then make the bwa index:

```bash
$ bwa index library.fa
```

And align the `test.cnt.fq` file produced by previous step to this library:

```bash
$ bwa aln library.fa test1.cnt.fq -n 0 > test1.sai  # run aln algorithm with no mismatch
$ bwa samse library.fa ./test1.sai ./test1.cnt.fq > test1.sam  # produce sam file
```

### 3. Recovery all bait-prey interaction pairs

Use tool `getedges`, for example:

```bash
$ ./getedges ./test1.cnt ./test1.sam -o test1.edges.tsv
```

More usage detail see:

```bash
./getedges -h
```

## About RLL-Y2H
More detail about RLL-Y2H please see the original paper:

```
Yang, Fang, et al. "Development and application of a recombination-based library versus library high-throughput yeast two-hybrid (RLL-Y2H) screening system." Nucleic acids research 46.3 (2018): e17-e17.
```

