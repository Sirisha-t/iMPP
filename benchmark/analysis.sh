#!/bin/bash
### BENCHMARK SCRIPTS ###

## Please update the paths with your local folder path 
SCRIPT="scripts"


## -i input reads (in fastq format)
## -p 0 or 1 ( 0= single, 1=paired)
## -g fgs/prodigal (type of genecaller)
## -t (comp/ref; comp=complete dataset, without reference. ref=subsampled dataset, with reference)
while getopts i:p:g:t: flag
do
    case "${flag}" in
        i) input=${OPTARG};;
        p) paired=${OPTARG};;
        g) genecaller=${OPTARG};;
        t) dstype=${OPTARG};;
    esac
done


if [ $genecaller = "fgs" ];
then
	bash $SCRIPT/runFGS_analysis.sh -i $input -p $paired -t $dstype
#elif [$genecaller == "mpd"];
#then
#	bash $SCRIPT/runMPD_analysis.sh -i $input -p $paired -t $dstype
#else
#	echo "Wrong genecaller parameter provided. Please use either 'fgs' or 'mpd'"
fi










