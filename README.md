# rust-fastp-clip

 - rust-fastp-clip
 - provide the reads in the fastq from novaseq, nextseq and clip the regions according to your choice.
 - a standalone application as well as a part of the rust-fastp package, faster than the fastp using tokio and async programming. 
 - no external dependencies, pure datastructure approach. 
 
 - build the release binary
 ```
 cargo build
 ```
 - run the compiled binary 

 ```
   ./rust-fastp-clip ./sample-files/samp2_L001_R1_001.fastq ./sample-files/samp2_L001_R2_001.fastq 0 10
 ```
 Gaurav Sablok
