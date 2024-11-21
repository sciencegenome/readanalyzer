use clap::Parser;

#[derive(Debug, Parser)]
#[clap(version)]

pub struct FastqArgs {
    /// please provide the reads R1 file path
    pub reads_1_arg: String,
    /// please provide the reads R2 file path
    pub reads_2_arg: String,
    /// please provide the clip region start
    pub clip_start: usize,
    /// please provide the clip region end
    pub clip_end: usize,
}
