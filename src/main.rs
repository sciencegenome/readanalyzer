mod args;
use args::FastqArgs;
use clap::Parser;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

/*
 *Author Gaurav Sablok
 *Universitat Potsdam
 *Date 2024-11-21

* implementing all the parts of the fastp in rustlang for the faster read user
* in nextseqseq, novaseq, and other high-throughput illumina platforms. These are
* available as separate rust-applications each of them and also a combined rust-fastp.
* This is the rust-fastp-clip, which will clip your reads at the specified regions.
*
* */

fn main() {
    let args = FastqArgs::parse();
    fastqclip(
        &args.reads_1_arg,
        &args.reads_2_arg,
        args.clip_start,
        args.clip_end,
    );

     fastq_quality_clip(
        &args.reads_1_arg,
        &args.reads_2_arg,
        args.clip_start,
        args.clip_end,
    );
}

fn fastqclip (fastq1: &str, fastq2: &str, clipregion: usize, cliplength: usize ) {

    #[derive(Debug, Clone)]
    struct FileFastqPre {
     header: String,
     sequence: String,
    }

    #[derive(Debug, Clone)]
    struct FileFastqPost {
      header: String,
      sequence: String,
      clippedregion: String
    }

    let file1 = File::open(fastq1).expect("file not present");
    let file2 = File::open(fastq2).expect("file not present");
    let fileread_1 = BufReader::new(&file1);
    let fileread_2 = BufReader::new(&file2);
    let mut fastq_1: Vec<FileFastqPre> = Vec::new();
    let mut fastq_2: Vec<FileFastqPre> = Vec::new();
    let mut fastq_1_header: Vec<String> = Vec::new();
    let mut fastq_2_header: Vec<String> = Vec::new();
    let mut fastq_1_sequence: Vec<String> = Vec::new();
    let mut fastq_2_sequence: Vec<String> = Vec::new();

    for i in fileread_1.lines(){
    let line = i.expect("line not present");
    if line.starts_with("@"){
        fastq_1_header.push(line);
    } else if line.starts_with("A") && !line.contains("E") || line.starts_with("T")
        && !line.contains("E") || line.starts_with("G") && !line.contains("E") ||
            line.starts_with("C") && !line.contains("E") || line.starts_with("N") && !line.contains("E")
    {
    fastq_1_sequence.push(line);
   }
}

    for i in fileread_2.lines(){
    let line = i.expect("line not present");
    if line.starts_with("@"){
        fastq_2_header.push(line);
    } else if line.starts_with("A") && !line.contains("E") || line.starts_with("T")
        && !line.contains("E") || line.starts_with("G") && !line.contains("E") ||
            line.starts_with("C") && !line.contains("E") || line.starts_with("N") && !line.contains("E")
    {
    fastq_2_sequence.push(line);
   }
}
   for i in 0..fastq_1_header.len() {
       fastq_1.push(FileFastqPre{
           header: fastq_1_header[i].clone(),
           sequence: fastq_1_sequence[i].clone()
       })
   }

    for i in 0..fastq_2_header.len() {
       fastq_2.push(FileFastqPre{
           header: fastq_2_header[i].clone(),
           sequence: fastq_2_sequence[i].clone(),
       })
   }

   let mut clipped_fastq_1: Vec<FileFastqPost> = Vec::new();
   let mut clipped_fastq_2: Vec<FileFastqPost> = Vec::new();

   for i in fastq_1.iter() {
       clipped_fastq_1.push(FileFastqPost{
           header: i.header.clone(),
           sequence: i.sequence.clone(),
           clippedregion: i.sequence[clipregion..cliplength].to_string(),
       })
   }

    for i in fastq_2.iter() {
       clipped_fastq_2.push(FileFastqPost{
           header: i.header.clone(),
           sequence: i.sequence.clone(),
           clippedregion: i.sequence[clipregion..cliplength].to_string(),
       })
   }

   let mut fastq_1_write = File::create("clipped-fastq-1.fastq").expect("file not present");
   let mut fastq_2_write = File::create("clipped-fastq-2.fastq").expect("file not found");

   for i in clipped_fastq_1.iter() {
       writeln!(fastq_1_write, "{}\n{}", i.header, i.clippedregion).expect("file not present");
   }
   for i in clipped_fastq_2.iter() {
       writeln!(fastq_2_write, "{}\n{}", i.header, i.clippedregion).expect("file not present");
   }

}

fn fastq_quality_clip (fastq1: &str, fastq2: &str, clipregion: usize, cliplength: usize ) {

    #[derive(Debug, Clone)]
    struct FileFastqPre {
     header: String,
     sequence: String,
     strand: String,
     quality:String,
    }

    #[derive(Debug, Clone)]
    struct FileFastqPost {
      header: String,
      sequence: String,
      clippedregion: String,
      strand:String,
      clippedquality:String,
    }


    let file1 = File::open(fastq1).expect("file not present");
    let file2 = File::open(fastq2).expect("file not present");
    let fileread_1 = BufReader::new(&file1);
    let fileread_2 = BufReader::new(&file2);
    let mut fastq_1: Vec<FileFastqPre> = Vec::new();
    let mut fastq_2: Vec<FileFastqPre> = Vec::new();
    let mut fastq_1_header: Vec<String> = Vec::new();
    let mut fastq_2_header: Vec<String> = Vec::new();
    let mut fastq_1_sequence: Vec<String> = Vec::new();
    let mut fastq_2_sequence: Vec<String> = Vec::new();
    let mut fastq_1_quality: Vec<String> = Vec::new();
    let mut fastq_2_quality: Vec<String> = Vec::new();
    let mut fastq_1_strand: Vec<String> = Vec::new();
    let mut fastq_2_strand: Vec<String> = Vec::new();

    for i in fileread_1.lines(){
    let line = i.expect("line not present");
    if line.starts_with("@"){
        fastq_1_header.push(line);
    } else if line.starts_with("A") && !line.contains("E") || line.starts_with("T")
        && !line.contains("E") || line.starts_with("G") && !line.contains("E") ||
            line.starts_with("C") && !line.contains("E") || line.starts_with("N") && !line.contains("E")
    {
    fastq_1_sequence.push(line);
    } else if line.starts_with("+") || line.starts_with("-"){
        fastq_1_strand.push(line);
    } else if line.contains("E"){
        fastq_1_quality.push(line)
    }
}

    for i in fileread_2.lines(){
    let line = i.expect("line not present");
    if line.starts_with("@"){
        fastq_2_header.push(line);
    } else if line.starts_with("A") && !line.contains("E") || line.starts_with("T")
        && !line.contains("E") || line.starts_with("G") && !line.contains("E") ||
            line.starts_with("C") && !line.contains("E") || line.starts_with("N") && !line.contains("E")
    {
      fastq_2_sequence.push(line);
   } else if line.starts_with("+") || line.starts_with("-"){
        fastq_2_strand.push(line);
    } else if line.contains("E"){
        fastq_2_quality.push(line)
    }
}
   for i in 0..fastq_1_header.len() {
       fastq_1.push(FileFastqPre{
           header: fastq_1_header[i].clone(),
           sequence: fastq_1_sequence[i].clone(),
           strand: fastq_1_strand[i].clone(),
           quality: fastq_1_quality[i].clone(),
       })
   }

    for i in 0..fastq_2_header.len() {
       fastq_2.push(FileFastqPre{
           header: fastq_2_header[i].clone(),
           sequence: fastq_2_sequence[i].clone(),
           strand: fastq_2_strand[i].clone(),
           quality: fastq_2_quality[i].clone(),
       })
   }

   let mut clipped_fastq_1: Vec<FileFastqPost> = Vec::new();
   let mut clipped_fastq_2: Vec<FileFastqPost> = Vec::new();

   for i in fastq_1.iter_mut() {
       clipped_fastq_1.push(FileFastqPost{
           header: i.header.clone(),
           sequence: i.sequence.clone(),
           clippedregion: i.sequence[clipregion..cliplength].to_string(),
           strand: i.strand.clone(),
           clippedquality: i.quality[clipregion..cliplength].to_string(),
       })
   }

      for i in fastq_2.iter() {
       clipped_fastq_2.push(FileFastqPost{
           header: i.header.clone(),
           sequence: i.sequence.clone(),
           clippedregion: i.sequence[clipregion..cliplength].to_string(),
           strand:i.strand.clone(),
           clippedquality: i.quality[clipregion..cliplength].to_string(),
       })
   }

   let mut fastq_1_write = File::create("clipped-fastqual-1.fastq")
   .expect("file not present");
   let mut fastq_2_write = File::create("clipped-fastqual-2.fastq")
   .expect("file not found");

   for i in clipped_fastq_1.iter() {
       writeln!(fastq_1_write, "{}\n{}\n{}\n{}", i.header, i.clippedregion, i.strand, i.clippedquality)
       .expect("file not present");
   }
   for j in clipped_fastq_2.iter() {
       writeln!(fastq_2_write, "{}\n{}\n{}\n{}", j.header, j.clippedregion, j.strand, j.clippedquality)
       .expect("file not present");
   }

}
