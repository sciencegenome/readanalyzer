use crate::filestruct::{
    Fastq1, Fastq1All, Fastq2, Fastq2All, SyncF1, SyncF1All, SyncF2, SyncF2All,
};
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

/*
Author Gaurav Sablok
SLB Potsdam
Date 2025-2-25

*/

pub fn fastqcompare(fastq1: &str, fastq2: &str) -> Result<String, Box<dyn Error>> {
    let open_f1 = File::open(&fastq1).unwrap();
    let open_f2 = File::open(&fastq2).unwrap();
    let read_f1 = BufReader::new(&open_f1);
    let read_f2 = BufReader::new(&open_f2);
    let mut header_f1 = Vec::new();
    let mut sequence_f1 = Vec::new();
    let mut header_f2 = Vec::new();
    let mut sequence_f2 = Vec::new();
    for i in read_f1.lines() {
        let line = i.expect("line not present");
        if line.starts_with("@") {
            header_f1.push(line);
        } else if line.starts_with("A") && !line.contains("E")
            || line.starts_with("T") && !line.contains("E")
            || line.starts_with("G") && !line.contains("E")
            || line.starts_with("C") && !line.contains("E")
            || line.starts_with("N") && !line.contains("E")
        {
            sequence_f1.push(line);
        }
    }
    for j in read_f2.lines() {
        let line = j.expect("line not present");
        if line.starts_with("@") {
            header_f2.push(line);
        } else if line.starts_with("A") && !line.contains("E")
            || line.starts_with("T") && !line.contains("E")
            || line.starts_with("G") && !line.contains("E")
            || line.starts_with("C") && !line.contains("E")
            || line.starts_with("N") && !line.contains("E")
        {
            sequence_f2.push(line);
        }
    }

    let mut combined_f1 = Vec::new();
    let mut combined_f2 = Vec::new();

    for i in 0..header_f1.len() {
        combined_f1.push(Fastq1 {
            header: header_f1[i].clone(),
            sequence: sequence_f1[i].clone(),
        });
    }
    for j in 0..header_f2.len() {
        combined_f2.push(Fastq2 {
            header: header_f2[j].clone(),
            sequence: sequence_f2[j].clone(),
        });
    }

    let mut sync_f1: Vec<SyncF1> = Vec::new();
    let mut sync_f2: Vec<SyncF2> = Vec::new();

    for i in combined_f1.iter() {
        for j in combined_f2.iter() {
            let f1_header = i.header.split(" ").collect::<Vec<&str>>()[0];
            let f2_header = j.header.split(" ").collect::<Vec<&str>>()[0];
            if f1_header.to_string() == f2_header.to_string() {
                sync_f1.push(SyncF1 {
                    header: f1_header.to_string(),
                    sequence: i.sequence.clone(),
                });
                sync_f2.push(SyncF2 {
                    header: f2_header.to_string(),
                    sequence: j.sequence.clone(),
                });
            }
        }
    }
    let mut f1_write = File::create("fasta1.fasta").unwrap();
    let mut f2_write = File::create("fasta2.fasta").unwrap();
    for i in sync_f1.iter_mut() {
        writeln!(f1_write, "{}\n{}", i.header, i.sequence).expect("not written");
    }
    println!("The read file for the first has been written");
    for j in sync_f2.iter_mut() {
        writeln!(f2_write, "{}\n{}", j.header, j.sequence).expect("not written");
    }
    println!("the read file for the second has been written");
    println!(
        "The total number of the synced reads in the nextseq/noavseq lanes are:{}",
        sync_f1.len()
    );
    Ok("The reads have been synced".to_string())
}

pub fn fastqgenerate(fastq1: &str, fastq2: &str) -> Result<String, Box<dyn Error>> {
    let open_f1 = File::open(&fastq1).unwrap();
    let open_f2 = File::open(&fastq2).unwrap();
    let read_f1 = BufReader::new(&open_f1);
    let read_f2 = BufReader::new(&open_f2);
    let mut header_f1 = Vec::new();
    let mut sequence_f1 = Vec::new();
    let mut strand_f1 = Vec::new();
    let mut quality_f1 = Vec::new();
    let mut header_f2 = Vec::new();
    let mut sequence_f2 = Vec::new();
    let mut strand_f2 = Vec::new();
    let mut quality_f2 = Vec::new();

    for i in read_f1.lines() {
        let line = i.expect("line not present");
        if line.starts_with("@") {
            header_f1.push(line);
        } else if line.starts_with("A") && !line.contains("E")
            || line.starts_with("T") && !line.contains("E")
            || line.starts_with("G") && !line.contains("E")
            || line.starts_with("C") && !line.contains("E")
            || line.starts_with("N") && !line.contains("E")
        {
            sequence_f1.push(line);
        } else if line.starts_with("+") || line.starts_with("-") {
            strand_f1.push(line);
        } else if line.contains("E") {
            quality_f1.push(line);
        }
    }
    for j in read_f2.lines() {
        let line = j.expect("line not present");
        if line.starts_with("@") {
            header_f2.push(line);
        } else if line.starts_with("A") && !line.contains("E")
            || line.starts_with("T") && !line.contains("E")
            || line.starts_with("G") && !line.contains("E")
            || line.starts_with("C") && !line.contains("E")
            || line.starts_with("N") && !line.contains("E")
        {
            sequence_f2.push(line);
        } else if line.starts_with("+") || line.starts_with("-") {
            strand_f2.push(line);
        } else if line.contains("E") {
            quality_f2.push(line);
        }
    }

    let mut combined_f1: Vec<Fastq1All> = Vec::new();
    let mut combined_f2: Vec<Fastq2All> = Vec::new();

    for i in 0..header_f1.len() {
        combined_f1.push(Fastq1All {
            header: header_f1[i].clone(),
            sequence: sequence_f1[i].clone(),
            strand: strand_f1[i].clone(),
            quality: quality_f1[i].clone(),
        });
    }

    for j in 0..header_f2.len() {
        combined_f2.push(Fastq2All {
            header: header_f2[j].clone(),
            sequence: sequence_f2[j].clone(),
            strand: strand_f2[j].clone(),
            quality: quality_f2[j].clone(),
        });
    }

    let mut sync_f1: Vec<SyncF1All> = Vec::new();
    let mut sync_f2: Vec<SyncF2All> = Vec::new();

    for i in combined_f1.iter() {
        for j in combined_f2.iter() {
            let f1_header = i.header.split(" ").collect::<Vec<&str>>()[0];
            let f2_header = j.header.split(" ").collect::<Vec<&str>>()[0];
            if f1_header.to_string() == f2_header.to_string() {
                sync_f1.push(SyncF1All {
                    header: f1_header.to_string(),
                    sequence: i.sequence.clone(),
                    strand: i.strand.clone(),
                    quality: i.quality.clone(),
                });
                sync_f2.push(SyncF2All {
                    header: f2_header.to_string(),
                    sequence: j.sequence.clone(),
                    strand: j.strand.clone(),
                    quality: j.quality.clone(),
                });
            }
        }
    }
    let mut f1_write = File::create("r1_analyzed.fastq").unwrap();
    let mut f2_write = File::create("r2_analyzed.fastq").unwrap();
    for i in sync_f1.iter_mut() {
        writeln!(
            f1_write,
            "{}\n{}\n{}\n{}",
            i.header, i.sequence, i.strand, i.quality
        )
        .expect("not written");
    }
    println!("The read file for the first has been written");
    for j in sync_f2.iter_mut() {
        writeln!(
            f2_write,
            "{}\n{}\n{}\n{}",
            j.header, j.sequence, j.strand, j.quality
        )
        .expect("not written");
    }
    println!("the read file for the second has been written");
    println!(
        "The total number of the synced reads in the nextseq/noavseq lanes are:{}",
        sync_f1.len()
    );
    Ok("The results have been synced".to_string())
}
