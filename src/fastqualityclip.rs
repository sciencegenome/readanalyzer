use crate::filestruct::FileFastqPostAdd;
use crate::filestruct::FileFastqPreAdd;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
/*
Author Gaurav Sablok
SLB Potsdam
Date 2025-2-25

*/

pub fn fastqualityclip(
    fastq1: &str,
    fastq2: &str,
    clipregion: usize,
    cliplength: usize,
) -> Result<String, Box<dyn Error>> {
    let file1 = File::open(fastq1).expect("file not present");
    let file2 = File::open(fastq2).expect("file not present");
    let fileread_1 = BufReader::new(&file1);
    let fileread_2 = BufReader::new(&file2);
    let mut fastq_1: Vec<FileFastqPreAdd> = Vec::new();
    let mut fastq_2: Vec<FileFastqPreAdd> = Vec::new();
    let mut fastq_1_header: Vec<String> = Vec::new();
    let mut fastq_2_header: Vec<String> = Vec::new();
    let mut fastq_1_sequence: Vec<String> = Vec::new();
    let mut fastq_2_sequence: Vec<String> = Vec::new();
    let mut fastq_1_quality: Vec<String> = Vec::new();
    let mut fastq_2_quality: Vec<String> = Vec::new();
    let mut fastq_1_strand: Vec<String> = Vec::new();
    let mut fastq_2_strand: Vec<String> = Vec::new();

    for i in fileread_1.lines() {
        let line = i.expect("line not present");
        if line.starts_with("@") {
            fastq_1_header.push(line);
        } else if line.starts_with("A") && !line.contains("E")
            || line.starts_with("T") && !line.contains("E")
            || line.starts_with("G") && !line.contains("E")
            || line.starts_with("C") && !line.contains("E")
            || line.starts_with("N") && !line.contains("E")
        {
            fastq_1_sequence.push(line);
        } else if line.starts_with("+") || line.starts_with("-") {
            fastq_1_strand.push(line);
        } else if line.contains("E") {
            fastq_1_quality.push(line)
        }
    }

    for i in fileread_2.lines() {
        let line = i.expect("line not present");
        if line.starts_with("@") {
            fastq_2_header.push(line);
        } else if line.starts_with("A") && !line.contains("E")
            || line.starts_with("T") && !line.contains("E")
            || line.starts_with("G") && !line.contains("E")
            || line.starts_with("C") && !line.contains("E")
            || line.starts_with("N") && !line.contains("E")
        {
            fastq_2_sequence.push(line);
        } else if line.starts_with("+") || line.starts_with("-") {
            fastq_2_strand.push(line);
        } else if line.contains("E") {
            fastq_2_quality.push(line)
        }
    }
    for i in 0..fastq_1_header.len() {
        fastq_1.push(FileFastqPreAdd {
            header: fastq_1_header[i].clone(),
            sequence: fastq_1_sequence[i].clone(),
            strand: fastq_1_strand[i].clone(),
            quality: fastq_1_quality[i].clone(),
        })
    }

    for i in 0..fastq_2_header.len() {
        fastq_2.push(FileFastqPreAdd {
            header: fastq_2_header[i].clone(),
            sequence: fastq_2_sequence[i].clone(),
            strand: fastq_2_strand[i].clone(),
            quality: fastq_2_quality[i].clone(),
        })
    }

    let mut clipped_fastq_1: Vec<FileFastqPostAdd> = Vec::new();
    let mut clipped_fastq_2: Vec<FileFastqPostAdd> = Vec::new();

    for i in fastq_1.iter_mut() {
        clipped_fastq_1.push(FileFastqPostAdd {
            header: i.header.clone(),
            sequence: i.sequence.clone(),
            clippedregion: i.sequence[clipregion..cliplength].to_string(),
            strand: i.strand.clone(),
            clippedquality: i.quality[clipregion..cliplength].to_string(),
        })
    }

    for i in fastq_2.iter() {
        clipped_fastq_2.push(FileFastqPostAdd {
            header: i.header.clone(),
            sequence: i.sequence.clone(),
            clippedregion: i.sequence[clipregion..cliplength].to_string(),
            strand: i.strand.clone(),
            clippedquality: i.quality[clipregion..cliplength].to_string(),
        })
    }

    let mut fastq_1_write = File::create("clipped-fastqual-1.fastq").expect("file not present");
    let mut fastq_2_write = File::create("clipped-fastqual-2.fastq").expect("file not found");

    for i in clipped_fastq_1.iter() {
        writeln!(
            fastq_1_write,
            "{}\n{}\n{}\n{}",
            i.header, i.clippedregion, i.strand, i.clippedquality
        )
        .expect("file not present");
    }
    for j in clipped_fastq_2.iter() {
        writeln!(
            fastq_2_write,
            "{}\n{}\n{}\n{}",
            j.header, j.clippedregion, j.strand, j.clippedquality
        )
        .expect("file not present");
    }

    Ok("The reads have been clipped".to_string())
}
