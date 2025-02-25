mod args;
mod clipall;
mod fastqualityclip;
mod filestruct;
use crate::clipall::
use crate::args::CommandParse;
use crate::args::Commands;
use crate::clipall::r1_process;
use crate::clipall::r2_process;
use crate::fastqualityclip::fastqualityclip;

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
        Commands::ClipFastqArgs {
            reads_1_arg,
            reads_2_arg,
            clip_start,
            clip_end,
        } => {
            let command = fastqualityclip(reads_1_arg, reads_2_arg, clip_start, clip_end).unwrap();

            println!("The commad has been completed:{:?}", command);
        }
        Commands::QualityClipArgs {
            reads_1_arg,
            reads_2_arg,
            quality_score,
        } => {
            let command1 = r1_process(reads_1_arg, quality_score ).unwrap();

           let command2 = r2_process(reads_1_arg, quality_score ).unwrap();
            println!("The command has been completed:{:?}", command1);
            println!("The command has been completed:{:?}", command2);
        }
    }
}
