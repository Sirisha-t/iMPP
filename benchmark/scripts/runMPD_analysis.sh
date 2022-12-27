## Please update the paths with your local folder path 
REF="test/ref.fna"
OUTDIR="output"
SCRIPT="scripts"
NF_DIR=".."
MPD_DIR="lib/prodigal"
SGA_DIR="lib"
SPADES_DIR="lib/spades/bin"
BWA_DIR="lib/bwa"
DIAMOND_DIR="lib/diamond"
PLASS_DIR="lib/plass/bin"
SPADES_OUT="output/spades"
CONFIG="config/config.txt"
threads=16

## read input arguments
while getopts i:p:t: flag
do
    case "${flag}" in
        i) input=${OPTARG};;
        p) paired=${OPTARG};;
        t) type=${OPTARG};;
    esac
done

## check if output dir exists, if not create one
if [ ! -d $OUTDIR ] 
then
    mkdir $OUTDIR
fi
##################################
## iMPP (MPD) run 
##################################
if [[$paired -eq 1]]
then
    $NF_DIR/nextflow $NF_DIR/main.nf --interleave $input --genecaller prodigal --outdir $OUTDIR -profile base,docker -with-report $OUTDIR/report.html -with-trace $OUTDIR/trace.txt -with-timeline $OUTDIR/timeline.html -resume &> $OUTDIR/impp_run.log
else
     $NF_DIR/nextflow $NF_DIR/main.nf --single $input --genecaller prodigal --outdir $OUTDIR -profile base,docker -with-report $OUTDIR/report.html -with-trace $OUTDIR/trace.txt -with-timeline $OUTDIR/timeline.html -resume &> $OUTDIR/impp_run.log
fi

###################################
##  MPD on reads  
###################################
## Command to convert fastq input to fasta format
## Note: MPD only accepts fasta input. (in case of fastq input, need to convert to fasta format)
sed -n '1~4s/^@/>/p;2~4p' $input > $OUTDIR/reads.fasta 
## run MPD on reads
$MPD_DIR/prodigal -i $OUTDIR/reads.fasta -d $OUTDIR/reads.ffn -f gff -o $OUTDIR/reads.gff -p meta -q -s $OUTDIR/reads.genes -a $OUTDIR/reads.faa


###################################
##  SGA run  
###################################
if [[$paired -eq 1]]
then
    $SGA_DIR/sga preprocess -o assembly.pp.fq -p 2 $input
else    
    $SGA_DIR/sga preprocess -o assembly.pp.fq $input
fi
$SGA_DIR/sga index --no-reverse -t $threads
$SGA_DIR/sga correct -t $threads -k 15 --learn -m 10 assembly.pp.fq -o assembly.ec.fq
$SGA_DIR/sga index -t $threads -a ropebwt assembly.ec.fq
$SGA_DIR/sga rmdup -t $threads assembly.ec.fq
$SGA_DIR/sga overlap -t $threads -m 20 assembly.ec.rmdup.fa
gunzip -f assembly.ec.rmdup.asqg.gz
mv assembly.ec.rmdup.asqg $OUTDIR
rm -f assembly.*

##Getting sga contigs, MPD predictions and bwa mapping to reads
$SGA_DIR/sga assemble -m 20 --transitive-reduction -g 0 -o $OUTDIR/sga $OUTDIR/assembly.ec.rmdup.asqg
$MPD_DIR/prodigal -i $OUTDIR/sga-contigs.fa -d $OUTDIR/contigs.ffn -f gff -o $OUTDIR/contigs.gff -q -s $OUTDIR/contigs.genes -a $OUTDIR/contigs.faa

##running bwa against sga
$BWA_DIR/bwa index -p $OUTDIR/sga.index $OUTDIR/contigs.ffn
$BWA_DIR/bwa mem $OUTDIR/sga.index $OUTDIR/reads.fastq -a -t 8 -k 19 -o $OUTDIR/sga.contig.sam


##################################
## SPAdes run 
##################################
if [[$paired -eq 1]]
then
    $SPADES_DIR/spades.py -t $threads -m 500 -o $OUTDIR --12 $input --meta --only-assembler -k auto 
else    
    $SPADES_DIR/spades.py -t $threads -m 500 -o $OUTDIR --single $input --only-assembler -k auto
fi
$SCRIPT/renameDB.py $OUTDIR/assembly_graph.fastg

##Getting MPD predictions of spades contigs, and bwa mapping to reads
$MPD_DIR/prodigal -i $OUTDIR/spades/contigs.fasta -d $OUTDIR/spades_contigs.ffn -f gff -o $OUTDIR/spades_contigs.gff -q -s $OUTDIR/spades_contigs.genes -a $OUTDIR/spades_contigs.faa
##running bwa against spades
$BWA_DIR/bwa index -p $OUTDIR/spa.index $OUTDIR/spades_contigs.ffn
$BWA_DIR/bwa mem $OUTDIR/spa.index $OUTDIR/reads.fastq -a -t $threads -k 19 -T 30 -o $dir/spades.contig.sam


####################################
## Ground truth generation ((ONLY FOR SUBSAMPLED))
####################################
if [[ $type == "ref"]]
then
## MPD on reference
$MPD_DIR/prodigal -i $REF -d $OUTDIR/ref.ffn -f gff -o $OUTDIR/ref.gff -q -s $OUTDIR/ref.genes -a $OUTDIR/ref.faa
## Mapping reads to FGS predicted reference
$BWA_DIR/bwa index -p $OUTDIR/ref $REF
$BWA_DIR/bwa mem $OUTDIR/ref $OUTDIR/reads.fastq -t $threads -k 19 -T 30 -o $OUTDIR/ref.sam
## Compute precision, recall and F1-score
## for iMPP(MPD), SGA+MPD, SPAdes+MPD and MPD
##NOTE: Update the config file provided with your local path folders
$SCRIPT/bin/eval $CONFIG &> $OUTDIR/precision_recall.log
fi

############################################
## Compute peptide assembly statistics (for subsampled or simulated datasets)
## for iMPP(MPD), MPD+PLASS and PLASS
###########################################
## Run PLASS on all input reads
$PLASS_DIR/plass assemble --threads $threads --min-length 20 --num-iterations 12 $OUTDIR/reads.fasta $OUTDIR/assembled_proteins.reads.faa $OUTDIR/plass.tmp
## Run PLASS on MPD predicted reads
$PLASS_DIR/plass assemble --threads $threads --min-length 20 --num-iterations 12 $OUTDIR/reads.ffn $OUTDIR/assembled_proteins.mpd.faa $OUTDIR/mpd.tmp
##delete tmp directories after plass run (take up too much space)
rm -rf $OUTDIR/plass.tmp $OUTDIR/mpd.tmp

##Filter out plass contigs <60aa
$BIN/filterPlassContig.py $OUTDIR/assembled_proteins.reads.faa
$BIN/filterPlassContig.py $OUTDIR/assembled_proteins.mpd.faa

##Run diamond makedb on all contigs
$DIAMOND_DIR/diamond makedb --in $OUTDIR/assembled_proteins.impp.re.60.faa -d $OUTDIR/impp_contig.index
$DIAMOND_DIR/diamond makedb --in $OUTDIR/assembled_proteins.reads.re.60.faa -d $OUTDIR/read_contig.index
$DIAMOND_DIR/diamond makedb --in $OUTDIR/assembled_proteins.mpd.re.60.faa -d $OUTDIR/mpd_contig.index

###Run diamond blastx on contigs against reads
$DIAMOND_DIR/diamond blastx -p $threads -d $OUTDIR/impp_contig.index -q $OUTDIR/reads.impp.fasta -o $OUTDIR/plass.read.impp.dmd
$DIAMOND_DIR/diamond blastx -p $threads -d $OUTDIR/read_contig.index -q $OUTDIR/reads.fasta -o $OUTDIR/plass.read.read.dmd
$DIAMOND_DIR/diamond blastx -p $threads -d $OUTDIR/mpd_contig.index -q $OUTDIR/reads.ffn -o $OUTDIR/plass.read.mpd.dmd


if [[ $type == "ref"]]
then
##Run diamond makedb on ground truth reference proteins (for subsampled and simulated)
$DIAMOND_DIR/diamond makedb --in $REFDIR/ref.faa -d $OUTDIR/ref_protein.index
elif [[ type == "comp"]]
$DIAMOND_DIR/diamond makedb --in $REFDIR/uniprot.faa -d $OUTDIR/ref_protein.index
else
echo "Incorrect type provided. Please select either 'complete' or 'ref' (for subsampled and simulated) \n\n";
fi


##Run diamond blastp on reference against contigs
$DIAMOND_DIR/diamond blastp -p $threads -d $OUTDIR/ref_protein.index -q $OUTDIR/assembled_proteins.impp.re.60.faa -o $dir/plass.contig.impp.dmd
$DIAMOND_DIR/diamond blastp -p $threads -d $OUTDIR/ref_protein.index -q $OUTDIR/assembled_proteins.reads.re.60.faa -o $dir/plass.contig.read.dmd
$DIAMOND_DIR/diamond blastp -p $threads -d $OUTDIR/ref_protein.index -q $OUTDIR/assembled_proteins.mpd.re.60.faa -o $dir/plass.contig.fgs.dmd

##Run diamond blastx of ref against impp reads, all reads, mpd reads
$DIAMOND_DIR/diamond blastx -p $threads -d $OUTDIR/ref_protein.index -q $OUTDIR/reads.impp.fasta -o $OUTDIR/plass.ref.impp.dmd
$DIAMOND_DIR/diamond blastx -p $threads -d $OUTDIR/ref_protein.index -q $OUTDIR/reads.fasta -o $OUTDIR/plass.ref.read.dmd
$DIAMOND_DIR/diamond blastx -p $threads -d $OUTDIR/ref_protein.index -q $OUTDIR/reads.ffn -o $OUTDIR/plass.ref.mpd.dmd


##Run script to get contig-level and read-level specificities and peptide statistics
#1. iMPP (MPD)
$SCRIPT/bin/pepstats $OUTDIR/plass.ref.impp.dmd $OUTDIR/plass.read.impp.dmd $OUTDIR/plass.contig.impp.dmd $OUTDIR/assembled_proteins.impp.re.60.faa $OUTDIR/ref.faa 60
#2. PLASS
$SCRIPT/bin/pepstats $OUTDIR/plass.ref.read.dmd $OUTDIR/plass.read.read.dmd $OUTDIR/plass.contig.read.dmd $OUTDIR/assembled_proteins.reads.re.60.faa $OUTDIR/ref.faa 60
#3. MPD+PLASS
$SCRIPT/bin/pepstats $OUTDIR/plass.ref.fgs.dmd $OUTDIR/plass.read.fgs.dmd $OUTDIR/plass.contig.fgs.dmd $OUTDIR/assembled_proteins.fgs.re.60.faa $OUTDIR/ref.faa 60

##Concatenate all peptide results into single file
echo "--------------iMPP(MPD) peptide stats-------------" > $OUTDIR/peptideStats_benchmark_results.txt
cat $OUTDIR/plass.read.impp.seqlen.60.count >> $OUTDIR/peptideStats_benchmark_results.txt
echo "--------------------------------------------------" >> $OUTDIR/peptideStats_benchmark_results.txt
echo "--------------MPD+PLASS peptide stats-------------" >> $OUTDIR/peptideStats_benchmark_results.txt
cat $OUTDIR/plass.read.mpd.seqlen.60.count >> $OUTDIR/peptideStats_benchmark_results.txt
echo "--------------------------------------------------" >> $OUTDIR/peptideStats_benchmark_results.txt
echo "--------------PLASS peptide stats-----------------" >> $OUTDIR/peptideStats_benchmark_results.txt
cat $OUTDIR/plass.read.read.seqlen.60.count >> $OUTDIR/peptideStats_benchmark_results.txt
echo "--------------------------------------------------" >> $OUTDIR/peptideStats_benchmark_results.txt









