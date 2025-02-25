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
