#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct FileFastqPreAdd {
    pub header: String,
    pub sequence: String,
    pub strand: String,
    pub quality: String,
}

#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct FileFastqPostAdd {
    pub header: String,
    pub sequence: String,
    pub clippedregion: String,
    pub strand: String,
    pub clippedquality: String,
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct FileFastqPre {
    pub header: String,
    pub sequence: String,
    pub strand: String,
    pub quality: String,
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct FastqWindow1 {
    pub header: String,
    pub sequence: String,
    pub strand: String,
    pub quality: String,
    pub start: usize,
    pub end: usize,
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct FastqWindow2 {
    pub header: String,
    pub sequence: String,
    pub strand: String,
    pub quality: String,
    pub start: usize,
    pub end: usize,
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct Fastq1 {
    pub header: String,
    pub sequence: String,
}
#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct Fastq2 {
    pub header: String,
    pub sequence: String,
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct SyncF1 {
    pub header: String,
    pub sequence: String,
}
#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct SyncF2 {
    pub header: String,
    pub sequence: String,
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct Fastq1All {
    pub header: String,
    pub sequence: String,
    pub strand: String,
    pub quality: String,
}
#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct Fastq2All {
    pub header: String,
    pub sequence: String,
    pub strand: String,
    pub quality: String,
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct SyncF1All {
    pub header: String,
    pub sequence: String,
    pub strand: String,
    pub quality: String,
}
#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct SyncF2All {
    pub header: String,
    pub sequence: String,
    pub strand: String,
    pub quality: String,
}
