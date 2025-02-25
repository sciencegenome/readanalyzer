use crate::filestruct::FileFastqPre;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

/*
Author Gaurav Sablok
SLB Potsdam
Date 2025-2-25

*/

pub fn fastq_quality_drop(
    fastq1: &str,
    fastq2: &str,
    quality: usize,
) -> Result<String, Box<dyn Error>> {
    let dropval: Vec<_> = vec![
        '(', ')', '*', '+', ',', '-', '.', '/', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
        ':', ';', '<', '=', '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K',
    ];
    let dropvec: Vec<_> = (7..43).collect::<Vec<_>>();
    let mut searched = 'a';
    for i in 0..dropval.len() {
        if dropvec[i] == quality {
            searched = dropval[i];
        }
    }

    let file1 = File::open(&fastq1).expect("file not present");
    let fileread_1 = BufReader::new(&file1);
    let mut fastq_1: Vec<FileFastqPre> = Vec::new();
    let mut fastq_1_header: Vec<String> = Vec::new();
    let mut fastq_1_sequence: Vec<String> = Vec::new();
    let mut fastq_1_strand: Vec<String> = Vec::new();
    let mut fastq_1_quality: Vec<String> = Vec::new();

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
        } else if line.starts_with("+") {
            fastq_1_strand.push(line)
        } else if line.contains("E") {
            fastq_1_quality.push(line)
        }
    }
    for i in 0..fastq_1_header.len() {
        fastq_1.push(FileFastqPre {
            header: fastq_1_header[i].clone(),
            sequence: fastq_1_sequence[i].clone(),
            strand: fastq_1_strand[i].clone(),
            quality: fastq_1_quality[i].clone(),
        })
    }

    let file2 = File::open(&fastq2).expect("file not present");
    let fileread_2 = BufReader::new(&file2);
    let mut fastq_2: Vec<FileFastqPre> = Vec::new();
    let mut fastq_2_header: Vec<String> = Vec::new();
    let mut fastq_2_sequence: Vec<String> = Vec::new();
    let mut fastq_2_strand: Vec<String> = Vec::new();
    let mut fastq_2_quality: Vec<String> = Vec::new();

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
        } else if line.starts_with("+") {
            fastq_2_strand.push(line)
        } else if line.contains("E") {
            fastq_2_quality.push(line)
        }
    }
    for i in 0..fastq_2_header.len() {
        fastq_2.push(FileFastqPre {
            header: fastq_2_header[i].clone(),
            sequence: fastq_2_sequence[i].clone(),
            strand: fastq_2_strand[i].clone(),
            quality: fastq_2_quality[i].clone(),
        })
    }

    let mut writecleaned_header_1: Vec<String> = Vec::new();
    let mut writecleaned_sequence_1: Vec<String> = Vec::new();
    let mut writecleaned_strand_1: Vec<String> = Vec::new();
    let mut writecleaned_quality_1: Vec<String> = Vec::new();
    for i in fastq_1.iter_mut() {
        let tupleseq_1 = i.sequence.clone().chars().collect::<Vec<_>>();
        let tuplequality_1 = i.quality.clone().chars().collect::<Vec<_>>();
        let mut tupleclean_1 = Vec::new();
        let mut tuplequalityclean_1 = Vec::new();
        let mut finaltupleclean_1 = Vec::new();
        let mut finaltuplequalityclean_1 = Vec::new();
        for j in 0..tupleseq_1.len() {
            if tuplequality_1[j] == searched {
                continue;
            } else {
                tupleclean_1.push(tupleseq_1[j].to_string());
                tuplequalityclean_1.push(tuplequality_1[j].to_string());
                finaltupleclean_1.push(tupleclean_1.join("").to_string());
                finaltuplequalityclean_1.push(tuplequalityclean_1.join("").to_string());
            }
        }
        writecleaned_header_1.push(i.header.to_string());
        writecleaned_sequence_1.push(tupleclean_1.join("").to_string());
        writecleaned_strand_1.push(i.strand.to_string());
        writecleaned_quality_1.push(tuplequalityclean_1.join("").to_string());
    }
    let mut fileopen_1 = File::create("single-quality-drop1.fastq").expect("file not present");
    for i in 0..writecleaned_header_1.len() {
        write!(
            fileopen_1,
            "{}\n{}\n{}\n{}\n",
            writecleaned_header_1[i],
            writecleaned_sequence_1[i],
            writecleaned_strand_1[i],
            writecleaned_quality_1[i]
        )
        .expect("file not found");
    }

    let mut writecleaned_header_2: Vec<String> = Vec::new();
    let mut writecleaned_sequence_2: Vec<String> = Vec::new();
    let mut writecleaned_strand_2: Vec<String> = Vec::new();
    let mut writecleaned_quality_2: Vec<String> = Vec::new();
    for i in fastq_2.iter_mut() {
        let tupleseq_2 = i.sequence.clone().chars().collect::<Vec<_>>();
        let tuplequality_2 = i.quality.clone().chars().collect::<Vec<_>>();
        let mut tupleclean_2 = Vec::new();
        let mut tuplequalityclean_2 = Vec::new();
        let mut finaltupleclean_2 = Vec::new();
        let mut finaltuplequalityclean_2 = Vec::new();
        for j in 0..tupleseq_2.len() {
            if tuplequality_2[j] == searched {
                continue;
            } else {
                tupleclean_2.push(tupleseq_2[j].to_string());
                tuplequalityclean_2.push(tuplequality_2[j].to_string());
                finaltupleclean_2.push(tupleclean_2.join("").to_string());
                finaltuplequalityclean_2.push(tuplequalityclean_2.join("").to_string());
            }
        }
        writecleaned_header_2.push(i.header.to_string());
        writecleaned_sequence_2.push(tupleclean_2.join("").to_string());
        writecleaned_strand_2.push(i.strand.to_string());
        writecleaned_quality_2.push(tuplequalityclean_2.join("").to_string());
    }
    let mut fileopen_2 = File::create("single-quality-drop2.fastq").expect("file not present");
    for i in 0..writecleaned_header_2.len() {
        write!(
            fileopen_2,
            "{}\n{}\n{}\n{}\n",
            writecleaned_header_2[i],
            writecleaned_sequence_2[i],
            writecleaned_strand_2[i],
            writecleaned_quality_2[i]
        )
        .expect("file not found");
    }
    Ok("The reads have been qualit dropped".to_string())
}

pub fn fastqclipremain(
    fastq1: &str,
    fastq2: &str,
    quality: usize,
) -> Result<String, Box<dyn Error>> {
    let dropval: Vec<_> = vec![
        '(', ')', '*', '+', ',', '-', '.', '/', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
        ':', ';', '<', '=', '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K',
    ];
    let dropvec: Vec<_> = (7..43).collect::<Vec<_>>();
    let mut searched = 'a';
    for i in 0..dropval.len() {
        if dropvec[i] == quality {
            searched = dropval[i];
        }
    }

    let file1 = File::open(&fastq1).expect("file not present");
    let fileread_1 = BufReader::new(&file1);
    let mut fastq_1: Vec<FileFastqPre> = Vec::new();
    let mut fastq_1_header: Vec<String> = Vec::new();
    let mut fastq_1_sequence: Vec<String> = Vec::new();
    let mut fastq_1_strand: Vec<String> = Vec::new();
    let mut fastq_1_quality: Vec<String> = Vec::new();

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
        } else if line.starts_with("+") {
            fastq_1_strand.push(line)
        } else if line.contains("E") {
            fastq_1_quality.push(line)
        }
    }
    for i in 0..fastq_1_header.len() {
        fastq_1.push(FileFastqPre {
            header: fastq_1_header[i].clone(),
            sequence: fastq_1_sequence[i].clone(),
            strand: fastq_1_strand[i].clone(),
            quality: fastq_1_quality[i].clone(),
        })
    }

    let file2 = File::open(&fastq2).expect("file not present");
    let fileread_2 = BufReader::new(&file2);
    let mut fastq_2: Vec<FileFastqPre> = Vec::new();
    let mut fastq_2_header: Vec<String> = Vec::new();
    let mut fastq_2_sequence: Vec<String> = Vec::new();
    let mut fastq_2_strand: Vec<String> = Vec::new();
    let mut fastq_2_quality: Vec<String> = Vec::new();

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
        } else if line.starts_with("+") {
            fastq_2_strand.push(line)
        } else if line.contains("E") {
            fastq_2_quality.push(line)
        }
    }
    for i in 0..fastq_2_header.len() {
        fastq_2.push(FileFastqPre {
            header: fastq_2_header[i].clone(),
            sequence: fastq_2_sequence[i].clone(),
            strand: fastq_2_strand[i].clone(),
            quality: fastq_2_quality[i].clone(),
        })
    }

    let mut writecleaned_header_1: Vec<String> = Vec::new();
    let mut writecleaned_sequence_1: Vec<String> = Vec::new();
    let mut writecleaned_strand_1: Vec<String> = Vec::new();
    let mut writecleaned_quality_1: Vec<String> = Vec::new();
    for i in fastq_1.iter_mut() {
        let tupleseq_1 = i.sequence.clone().chars().collect::<Vec<_>>();
        let tuplequality_1 = i.quality.clone().chars().collect::<Vec<_>>();
        let mut tupleclean_1 = Vec::new();
        let mut tuplequalityclean_1 = Vec::new();
        let mut finaltupleclean_1 = Vec::new();
        let mut finaltuplequalityclean_1 = Vec::new();
        for j in 0..tupleseq_1.len() {
            for dropvalue in 0..dropval.len() {
                if tuplequality_1[j] == searched || tuplequality_1[j] != dropval[dropvalue] {
                    continue;
                } else {
                    tupleclean_1.push(tupleseq_1[j].to_string());
                    tuplequalityclean_1.push(tuplequality_1[j].to_string());
                    finaltupleclean_1.push(tupleclean_1.join("").to_string());
                    finaltuplequalityclean_1.push(tuplequalityclean_1.join("").to_string());
                }
            }
        }
        writecleaned_header_1.push(i.header.to_string());
        writecleaned_sequence_1.push(tupleclean_1.join("").to_string());
        writecleaned_strand_1.push(i.strand.to_string());
        writecleaned_quality_1.push(tuplequalityclean_1.join("").to_string());
    }
    let mut fileopen_1 = File::create("single-quality-all-drop1.fastq").expect("file not present");
    for i in 0..writecleaned_header_1.len() {
        write!(
            fileopen_1,
            "{}\n{}\n{}\n{}\n",
            writecleaned_header_1[i],
            writecleaned_sequence_1[i],
            writecleaned_strand_1[i],
            writecleaned_quality_1[i]
        )
        .expect("file not present");
    }

    let mut writecleaned_header_2: Vec<String> = Vec::new();
    let mut writecleaned_sequence_2: Vec<String> = Vec::new();
    let mut writecleaned_strand_2: Vec<String> = Vec::new();
    let mut writecleaned_quality_2: Vec<String> = Vec::new();
    for i in fastq_2.iter_mut() {
        let tupleseq_2 = i.sequence.clone().chars().collect::<Vec<_>>();
        let tuplequality_2 = i.quality.clone().chars().collect::<Vec<_>>();
        let mut tupleclean_2 = Vec::new();
        let mut tuplequalityclean_2 = Vec::new();
        let mut finaltupleclean_2 = Vec::new();
        let mut finaltuplequalityclean_2 = Vec::new();
        for j in 0..tupleseq_2.len() {
            for dropvalue in 0..dropval.len() {
                if tuplequality_2[j] == searched || tuplequality_2[j] != dropval[dropvalue] {
                    continue;
                } else {
                    tupleclean_2.push(tupleseq_2[j].to_string());
                    tuplequalityclean_2.push(tuplequality_2[j].to_string());
                    finaltupleclean_2.push(tupleclean_2.join("").to_string());
                    finaltuplequalityclean_2.push(tuplequalityclean_2.join("").to_string());
                }
            }
        }
        writecleaned_header_2.push(i.header.to_string());
        writecleaned_sequence_2.push(tupleclean_2.join("").to_string());
        writecleaned_strand_2.push(i.strand.to_string());
        writecleaned_quality_2.push(tuplequalityclean_2.join("").to_string());
    }
    let mut fileopen_2 = File::create("single-quality-all-drop2.fastq").expect("file not present");
    for i in 0..writecleaned_header_2.len() {
        write!(
            fileopen_2,
            "{}\n{}\n{}\n{}\n",
            writecleaned_header_2[i],
            writecleaned_sequence_2[i],
            writecleaned_strand_2[i],
            writecleaned_quality_2[i]
        )
        .expect("file not found");
    }

    Ok("The reads have been dropped".to_string())
}
