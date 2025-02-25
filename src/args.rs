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
    /// clipping specific regions
    ClipFastq {
        /// please provide the reads R1 file path
        reads_1: String,
        /// please provide the reads R2 file path
        reads_2: String,
        /// please provide the clip region start
        clip_start: usize,
        /// please provide the clip region end
        clip_end: usize,
    },
    /// drop specific quality
    QualityClip {
        /// please provide the reads R1 file path
        reads_1: String,
        /// please provide the reads R2 file path
        reads_2: String,
        /// please provide the quality value to be used as a threshold
        quality_score: usize,
    },
    /// drop specific score
    DropQuality {
        /// please provide the reads R1 file path
        reads_1: String,
        /// please provide the reads R2 file path
        reads_2: String,
        /// please provide the quality value to be used as a threshold
        quality_score: usize,
    },
    /// remove the quality and the adapter
    AdapterClipper {
        /// please provide the reads R1 file path
        reads_1: String,
        /// please provide the reads R2 file path
        reads_2: String,
        /// please provide the quality value to be used as a threshold
        quality_score: usize,
        /// please provide the adapter sequence
        adapter: String,
    },
    Sync {
        /// please provide the reads R1 file path
        reads_1: String,
        /// please provide the reads R2 file path
        reads_2: String,
    },
}
