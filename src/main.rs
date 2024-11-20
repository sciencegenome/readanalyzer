mod args;
use args::FastqArgs;
use clap::Parser;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

/*
 *Author Gaurav Sablok
 *Universitat Potsdam
 *Date 2024-11-20

* implementing all the parts of the fastp in rustlang for the faster read user
* in nextseqseq, novaseq, and other high-throughput illumina platforms. These are
* available as separate rust-applications each of them and also a combined rust-fastp.
* This is the rust-fastp-clip, which will clip your reads at the specified regions.
*
* */

fn main() {
    let args = FastqArgs::parse();
    fastqcompare(
        &args.reads_1_arg,
        &args.reads_2_arg
    );
 fastqgenerate(&args.reads_1_arg, &args.reads_2_arg);
}

