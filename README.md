==========================================================================

# iMPP: integrated Metagenomic Protein Predictor #

==========================================================================
### Overview ###

iMPP is a tool designed to predict genes and assemble peptides from short fragmentary reads.
iMPP is written in C++ and has been tested on a 64-bit Linux system.

The input for the software are FASTQ sequences, and the output comprises of four files:
* predicted genes (nucl) 
* predicted peptides (prot) 
* gene predictions in gff format
* assembled peptide sequences (prot).

The workflow for iMPP has been written in [Nextflow](https://www.nextflow.io/index.html), since it supports the use of [Docker](https://www.docker.com/) containers. 

## Instructions to run iMPP using Nextflow (with Docker enabled) ##

### Installation ###
__Install Java__

Make sure that Java is installed. (11 or higher) (Installation instructions can be found [here]((https://www.oracle.com/java/technologies/downloads/))
```
java --version
```

__Install Docker__

Follow [Docker Documentation](https://docs.docker.com/get-docker/) (docker CE is sufficient). Instructions to install Docker for Windows can be found [here](https://docs.docker.com/desktop/install/windows-install/) and for Mac [here](https://docs.docker.com/desktop/install/mac-install/).

Also follow the post-installation step to manage Docker as a non-root user ([here](https://docs.docker.com/engine/install/linux-postinstall/) for Linux), otherwise you will need to change the sudo option in nextflow docker config scope as described in the [nextflow documentation](https://www.nextflow.io/docs/latest/config.html#scope-docker). Non-privileged users can run Docker following the instrcutions provided [here](https://docs.docker.com/engine/security/rootless/).

```
To install Docker using CLI (Linux):
curl -fsSL https://get.docker.com | sh;
[Note: In case of permission denied issue, type: sudo chmod 666 /var/run/docker.sock]
```


__Clone the repository__
```
$git clone https://github.com/Sirisha-t/iMPP.git

Go to the project directory (/iMPP):

$cd iMPP/
```

__Install Nextflow__
```
curl -fsSL https://get.nextflow.io | bash

Make the binary executable on your system by running chmod +x nextflow.

Optionally, move the nextflow file to a directory accessible by your $PATH variable. 
```

__Run iMPP using Nextflow__
```
For single end:
/path/to/nextflow/nextflow run main.nf --single <input_fastq_file> --outdir <output_dir> -profile base,docker

For paired end:
/path/to/nextflow/nextflow run main.nf --forward <forward_fastq_file> --reverse <reverse_fastq_file> --outdir <output_dir> -profile base,docker


Example:
for single end: 
/path/to/nextflow/nextflow run main.nf --single data/reads.fq.gz --outdir output -profile base,docker

for paired end:
/path/to/nextflow/nextflow run main.nf --forward data/read1.fq.gz --reverse data/read2.fq.gz --outdir output -profile base,docker

```
* The input parameters can be modified based on the parameter options provided below:
```
USAGE: nextflow run main.nf --single/--forward and --reverse/--interleaved [other options] <fastq file/s> --outdir <output_directory> -profile base,docker
 Input data options:
   -profile               <string>      : [required] docker and base/test
   -resume				: [optional] can be set to resume workflow execution
   --single               <filename>    : fastq file with unpaired reads
   --forward              <filename>    : fastq file with forward paired-end reads
   --reverse              <filename>    : fastq file with reverse paired-end reads
   --interleaved          <filename>    : fastq file with interlaced forward and reverse paired-end reads
   --outdir               <dirname>     : [required] output directory name including the full path [default: output]
   --genecaller           <string>      : [optional] fgs or prodigal [default: fgs]
   --threads              <int>         : [optional] number of threads [default: 16]
   --maxlen               <int>         : [optional] maximum extension length for anchors [default: 150]
   -h/--help                            :  help message

  NOTE: Input FASTQ read file must be specified in either --interleaved, --forward & --reverse or --single format.
```
* The defult parameters used to run iMPP can be found here: [parameters.config](https://github.com/Sirisha-t/iMPP/blob/master/params/parameters.txt "parameters.txt")


* You can run a test with the following command:
```
 /path/to/nextflow/nextflow run main.nf -profile test,docker
 
 [This will automatically read the sample input file from the /data folder]
 
```
* To resume workflow execution from where it stopped, use the '-resume' command as shown below:
```
/path/to/nextflow/nextflow run main.nf -profile test,docker -resume
```

## Output ##

The final iMPP output should contain four files.
```
1. orfs.ffn : This file lists nucleotide sequences.
E.g.
>1_+
GTTGTTACCTCGTTACCTTTGGTCGAAAAAAAAAGCCCGCACTGTCAGGTGCGGGCTTTTTTCTGTGTTTCCTGTACGCGTCAGCCCGCACCGTTACCTG
TGGTAATGGTGATGGTGGTGGTAATGGTGGTGCTAATGCGTTTCATGGATGTTGTGTACTCTGTAATTTTTATCTGTCTGTGCGCTATGCCTATATTGGT
TAAAGTATTTAGTGACCTAAGTCAA
>2,17,14_90_191_+
TTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCC
TCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTAT
TTTTGCCGAACTTTTGACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTTGCCCAAATAAAACAT
GTCCTGCATGGCATTAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACGCTGCGCTGATTTGCCGTGGCGAGAAAATGTCGATCGCCATTATGGCCGGCG

2. orfs.faa : This file lists the amino acid sequences of the nucleotide sequneces in "[out_file_name].ffn".
E.g.
>17,51,49,46,41_2_136_+
VVTSLPLVEKKSPHCQVRAFFCVSCTRQPAPLPVVMVMVVVMVVLMRFMDVVYSVIFICLCAMPILVKVFSDLSQ
>2,17,14_90_191_+
LKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDALPNISDAERIFAELLTGLAAAQPGFPLAQLKTFVDQEFAQIKH
VLHGISLLGQCPDSINAALICRGEKMSIAIMAGVLEARGHNVTVIDPVEKLLAVGHYLESTVDIAESTRRIAASRIPADHMVLMAGFTAGNEKGELVVLG
RNGSDYSAAVLAACLRADCCEIWTDVDGVYTCDPRQVPDARLLKSMSYQEAMELSYFGAKVLHPRTITPIAQFQIPCLIKNTGNPQAPGTLIGASRDEDE

3. orfs.gff : This file lists the gene prediction results in gff format.
E.g.
##gff-version 3
2,17,14			FGS	CDS	90	191	.	+	2	ID=2,17,14_90_191_+;product=predicted protein
3,17,15			FGS	CDS	90	194	.	+	2	ID=3,17,15_90_194_+;product=predicted protein
14,51,49,46,40,59,58	FGS	CDS	2	136	.	+	1	ID=14,51,49,46,40,59,58_2_136_+;product=predicted protein
15,73,40,59,58		FGS	CDS	1	132	.	+	0	ID=15,73,40,59,58_1_132_+;product=predicted protein

4. assembled_proteins.faa : This file containts the assembled protein sequences.
E.g.
>0	1+98	3
LRSNLLKDFQEVIDDSKLKVVRNGYNGEILEVPAEKR
>1	1+98	3
ATNFPSIVDSELIELITDLLPTRCLIDTQVFDEEGFYRM
>1	99-98	3
HLFVTIKEVSDNPVLHPIKTLFIEDLCVDQAARGQKIGDQLYQFAVNYAREIGCYNLTLNVWN
>6	99-98	3
AWELMLKAYIINNNGEESIYFKDSKDRTISLENAVE
```
