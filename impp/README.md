==============================================

# iMPP: integrated Metagenomic Protein Predictor #

==============================================
### Description ###

iMPP is a tool designed to predict and assemble peptides from short fragmentary reads.
iMPP is written in C++ and has been tested on a 64-bit Linux system.
The input for the software are FASTQ sequences, and the output comprises of 4 files:
1. predicted genes (nucl), 2. predicted peptides (prot), 3. gene predictions in gff format, 4. assembled peptide sequences (prot).

==============================================
### Prerequisites ###

1. gcc compiler (version > 4.8.5)
2. boost-1.54.0 or newer
3. perl interpreter

==============================================
## Instructions to run using Nextflow with Docker ##

Nextflow and Docker would need to be installed to run the software. 

The below command can be copied to the terminal to install nextflow:
```
curl -fsSL https://get.nextflow.io | bash
```
You can then run the command:
```
nextflow run tsirisha/impp-nf --12 reads.fq -o output_dir

```
The input parameters can be modified based on the parameter options provided below:
```
	 --12  		 <filename>    :	 fastq file with interlaced forward and reverse paired-end reads
 	 --1             <filename>    :	 fastq file with forward paired-end reads
 	 --2	 	 <filename>    :	 fastq file with reverse paired-end reads
 	 --s	  	 <filename>    :	 fastq file with unpaired reads
 	 --assembler     <0 or 1>      :         type of assembler to run (0: sga, 1: spades)
	 --genecaller	 <fgs or mpd>  :         type of gene caller to run
 	 --outdir        <dirname>     :         output directory name including the full path
 	 --max-len       <int>         :         maximum extension length for anchors (default: 300) 
 	 --param-file    <filename>    :         parameter file
```

The docker image impp-nf will can also be pulled before running the above command:
```
docker pull tsirisha/impp-nf
```

Note: If you are unable to run the script, please check your Docker or Nextflow installation. 




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
The defult parameters used to run iMPP can be found here: [default_parameters](https://github.com/Sirisha-t/iMPP/blob/master/impp/params/parameters.txt "parameters.txt")


__Output__

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
==============================================
