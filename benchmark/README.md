### Instructions to run benchmarking on datasets

__All benchmarking scripts can be found in the /benchmark folder and its subfolders__

First step is to compile benchmarking scripts (/scripts/src folder) and install all 3rd party tools used for benchmarking (/lib folder). Please run the following setup script to install and compile and necessary tools/scripts.
```
./setup.sh
```

__Command line to run benchmarking on a test dataset is provided below__

For more details on the test dataset, please refer to Supplementary Methods and Results (SDS1, 2x coverage)

For testing on other datasets benchmarked in the manuscript, please download the reads from: https://cbb.ittc.ku.edu/iMPP.html

For complete datasets, please use the SRA accession numbers provided in the main manuscript to download the data first

Wrapper script (analysis.sh) command line arguments:
```
     -i input reads (in fastq format)
     -p 0 or 1 (0= single, 1=paired)
     -g fgs/prodigal (type of genecaller)
     -t (comp/ref; comp=complete dataset, without reference. ref=subsampled dataset, with reference)
```

__Run the script 'analysis.sh' with the following command line arguments to benchmark with FGS__

NOTE: Before running the following command, please ensure you update the path to all third-party tools and input/output directories in ~/benchmark/scripts/runFGS_analysis.sh (update to include your local path)

The test dataset and ground truth reference sequence are provided in the /test folder. You can update the folder to include other datasets (reads and ground truth) to benchmark. 
```
bash analysis.sh -i test/sds1.2x.fq -p 0 -g fgs -t ref
```
 

__Run the script 'analysis.sh' with the following command line arguments to benchmark with Prodigal__

NOTE: Before running the follwing command, please ensure the path/to/scripts are all correctly specified in ~/benchmark/scripts/runMPD_analysis.sh (update to include your local path )
```
bash analysis.sh -i test/sds1.2x.fq -p 0 -g mpd -t ref
```



