# rust-fastp-clip

 - rust-fastp-clip
 - provide the reads in the fastq from novaseq, nextseq and clip the regions according to your choice.
 - no external dependencies, pure datastructure approach.
 - please see the last commit message and if it says compiled binary then it is completed or else still in development version.
 
 - build the release binary
 ```
 cargo build
 ```
 ```
 ╭─gauravsablok@fedora ~/Downloads/rust-fastp-clip-main  
 ╰─➤  ./rust-fastp-clip -h
Usage: rust-fastp-clip <READS_1_ARG> <READS_2_ARG> <CLIP_START> <CLIP_END>

Arguments:
  <READS_1_ARG>  please provide the reads R1 file path
  <READS_2_ARG>  please provide the reads R2 file path
  <CLIP_START>   please provide the clip region start
  <CLIP_END>     please provide the clip region end

Options:
  -h, --help     Print help
  -V, --version  Print version
 - run the compiled binary 
 ```

 ```
   ./rust-fastp-clip ./sample-files/samp2_L001_R1_001.fastq ./sample-files/samp2_L001_R2_001.fastq 0 10
 ```
 Gaurav Sablok
