==========================================================================

# iMPP: integrated Metagenomic Protein Predictor #

==========================================================================
### Overview ###

iMPP is a tool designed to predict genes and assemble peptides from short fragmentary reads.
iMPP is written in C++ and has been tested on a 64-bit Linux system.

The input for the software are FASTQ sequences, and the output comprises of four files:
1. predicted genes (nucl), 2. predicted peptides (prot), 3. gene predictions in gff format, 4. assembled peptide sequences (prot).

The workflow for iMPP has been written in [Nextflow](https://www.nextflow.io/index.html), since it supports the use of [Docker](https://www.docker.com/) containers. 

### Prerequisites ###

1. gcc compiler (version > 4.8.5)
2. boost-1.54.0 or newer
3. perl interpreter
4. python (version >= 2.7)
5. Bash 3.2 (or later)


## Instructions to run iMPP using Nextflow (with Docker enabled) ##

### Installation ###
__Install Java__

Make sure that [Java](https://www.java.com/en/download/) is installed. (8 or higher)
```
java --version
```

__Install Nextflow__
```
curl -fsSL https://get.nextflow.io | bash

Add Nextflow binary to your user's PATH:
mv nextflow ~/bin/
OR system-wide installation:
sudo mv nextflow /usr/local/bin
```

__Install Docker__

Follow [Docker Documentation](https://docs.docker.com/get-docker/) (docker CE is sufficient). Also follow the post-installation step to manage Docker as a non-root user ([here](https://docs.docker.com/engine/install/linux-postinstall/) for Linux), otherwise you will need to change the sudo option in nextflow docker config scope as described in the nextflow documentation [here](https://www.nextflow.io/docs/latest/config.html#scope-docker). 

```
To get Docker using CLI:
curl -fsSL https://get.docker.com | sh;

Note: In case of permission denied issue, type: sudo chmod 666 /var/run/docker.sock
```

### Steps to run iMPP ###

Step 1: Clone the repository:
```
$git clone https://github.com/Sirisha-t/iMPP.git
```
Step 2: Go to the project directory (/iMPP)
```
$cd iMPP
```
Step 3: Run Nextflow script
```
For single end:
nextflow run main.nf --single <input_fastq_file> --outdir <output_dir> -profile base,docker

For paired end:
nextflow run main.nf --forward <forward_fastq_file> --reverse <reverse_fastq_file> --outdir <output_dir> -profile base,docker
OR
nextflow run main.nf --interleaves <interleaved_fastq_file> --outdir <output_dir> -profile base,docker


Eg:
for single end: 
nextflow run main.nf --single data/reads.fq.gz --outdir output -profile base, docker

for paired end:
nextflow run main.nf --forward data/read1.fq.gz --reverse data/read2.fq.gz --outdir output -profile base, docker
OR
nextflow run main.nf --interleaved data/reads.12.fq.gz --outdir output -profile base, docker

```
The input parameters can be modified based on the parameter options provided below:
```
USAGE: nextflow run main.nf --single/--forward and --reverse/--interleaved [other options] <fastq file/s> --outdir <output_directory> -profile base,docker
 Input data options:
   -profile               <string>      : [required] docker and base/test
   -resume				: [optional] can be set to resume workflow execution 
   --interleaved          <filename>    : fastq file with interlaced forward and reverse paired-end reads
   --forward              <filename>    : fastq file with forward paired-end reads
   --reverse              <filename>    : fastq file with reverse paired-end reads
   --single               <filename>    : fastq file with unpaired reads
   --outdir               <dirname>     : [required] output directory name including the full path [default: output]
   --genecaller           <string>      : [optional] fgs or prodigal [default: fgs]
   --threads              <int>         : [optional] number of threads [default: 16]
   --maxlen               <int>         : [optional] maximum extension length for anchors [default: 150]
   -h/--help                            :  help message

  NOTE: Input FASTQ read file must be specified in either --interleaved, --forward & --reverse or --single format.
```
*The defult parameters used to run iMPP can be found here: [parameters.config](https://github.com/Sirisha-t/iMPP/params/parameters.txt "parameters.txt")


*Note: If you are unable to run the script, please check your Docker or Nextflow installation. 

*You can run a test with the following command:
```
 nextflow run main.nf -profile test,docker
 
 [This will automatically read the sample input file from the /data folder]
 
```
*To resume workflow execution from where it stopped, use the '-resume' command as shown below:
```
nextflow run main.nf -profile test,docker -resume
```


## Instructions to install locally ##

To install iMPP in your local dir, please follow the steps below:

1. Clone the repository:
   $git clone https://github.com/Sirisha-t/iMPP.git
   
   This will generate the directory "iMPP".
   

2. Install 3rd party softwares and compile by running the "install.sh" script.
    $ bash install.sh


__Running the program:__

1.  The impp_run.pl perl wrapper is used to run the program. The command used is:

```
USAGE: ./impp_run.pl [options] <fastq_file/s> -o <output_directory>

Input data options:
	 --12  		 <filename> :	 fastq file with interlaced forward and reverse paired-end reads
 	 -1    		 <filename> :	 fastq file with forward paired-end reads
 	 -2	 	 <filename> :	 fastq file with reverse paired-end reads
 	 -s	  	 <filename> :	 fastq file with unpaired reads
 	 -a/--assembler  <0 or 1>   :    type of assembler to run (0: sga, 1: spades)
 	 -o/--outdir     <dirname>  :    output directory name including the full path
 	 -m/--max-len    <int>      :    maximum extension length for anchors (default: 300) 
 	 -p/--param-file <filename> :    parameter file
 	 -h/--help                  :    print help message

Note: The parameter file specifies the parameters used for running all third party programs. 
The 'parameters.txt' file in ~/params/ directory contains the deaulft set of parameters used for running iMPP. 
This file can be edited and provided as an optional parameter depending on user requirements. If no file is specified, 
preset parameters will be used.  
```

__Example__

An example simulated single-end Illumina reads file is provided in the ~/example/ directory.

iMPP can be run on this file with the following options:
```
$ perl impp_run.pl -s example/samplereads.fq -a 0 -o example -m 300 -p params/parameters.txt &> example/impp.example.run.log

(Note: The ~/example/samplereads.fq.gz is a zipped file. Please make sure you unzip it (cmd: gunzip samplereads.fq.gz) before running iMPP)  
```

In this example, iMPP will use the preset parameters from the params/parameters.txt file.

The defult parameters used to run iMPP can be found here: [parameters.config](https://github.com/Sirisha-t/iMPP/params/parameters.txt "parameters.txt")


## Output ##

The final iMPP output should contain four files.
```
1. orfs.ffn : This file lists nucleotide sequences.
E.g.
>17,51,49,46,41_2_136_+
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
