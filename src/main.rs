mod args;
mod clipall;
mod clipbase;
mod fastqualityclip;
mod filestruct;
mod qualityadapter;
mod sync;
use crate::args::CommandParse;
use crate::args::Commands;
use crate::clipall::r1_process;
use crate::clipall::r2_process;
use crate::clipbase::fastq_quality_drop;
use crate::clipbase::fastqclipremain;
use crate::fastqualityclip::fastqualityclip;
use crate::qualityadapter::r1adapter_process;
use crate::qualityadapter::r2adapter_process;
use crate::sync::fastqcompare;
use crate::sync::fastqgenerate;

use clap::Parser;

/*
Author Gaurav Sablok
SLB Potsdam
Date 2025-2-25

readanalyzer for the illumina nextseq, miseq, and novaseq.


*/

fn main() {
    let args = CommandParse::parse();
    match &args.command {
        Commands::ClipFastq {
            reads_1,
            reads_2,
            clip_start,
            clip_end,
        } => {
            let command = fastqualityclip(reads_1, reads_2, *clip_start, *clip_end).unwrap();

            println!("The commad has been completed:{:?}", command);
        }
        Commands::QualityClip {
            reads_1,
            reads_2,
            quality_score,
        } => {
            let command1 = r1_process(reads_1, *quality_score).unwrap();

            let command2 = r2_process(reads_2, *quality_score).unwrap();
            println!("The command has been completed:{:?}", command1);
            println!("The command has been completed:{:?}", command2);
        }
        Commands::DropQuality {
            reads_1,
            reads_2,
            quality_score,
        } => {
            let command1 = fastq_quality_drop(reads_1, reads_2, *quality_score).unwrap();
            let command2 = fastqclipremain(reads_1, reads_2, *quality_score).unwrap();
            println!(
                "The command have been completed:{:?}-{:?}",
                command1, command2
            );
        }
        Commands::AdapterClipper {
            reads_1,
            reads_2,
            quality_score,
            adapter,
        } => {
            let fncall1 = r1adapter_process(reads_1, *quality_score, adapter).unwrap();
            let fncall2 = r2adapter_process(reads_2, *quality_score, adapter).unwrap();
            println!("The rustlang clipall-R1 has been finished:{:?}", fncall1);
            println!("The rustlang clipall-R2 has been finished:{:?}", fncall2);
        }
        Commands::Sync { reads_1, reads_2 } => {
            let command1 = fastqcompare(reads_1, reads_2).unwrap();

            let command2 = fastqgenerate(reads_1, reads_2).unwrap();
            println!("The reads have been synced:{:?}-{:?}", command1, command2);
        }
    }
}
