use clap::{Parser, Subcommand};
#[derive(Debug, Parser)]
#[command(
    name = "readanalyzer",
    version = "1.0",
    about = "readanalyzer: analyze reads from illumina and long reads"
)]
pub struct CommandParse {
    /// subcommands for the specific actions
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// convert into fasta
    ClipFastqArgs {
        /// please provide the reads R1 file path
        reads_1_arg: String,
        /// please provide the reads R2 file path
        reads_2_arg: String,
        /// please provide the clip region start
        clip_start: usize,
        /// please provide the clip region end
        clip_end: usize,
    },
    QualityClipArgs {
        /// please provide the reads R1 file path
        reads_1_arg: String,
        /// please provide the reads R2 file path
        reads_2_arg: String,
        /// please provide the quality value to be used as a threshold
        quality_score: usize,
    },
}
