#bin/bash
#
# Script to install third-party tools and compile benchmarking scripts
#------------------------
# Pre-requisites needed
#------------------------
# cmake 
# zlib 
# gcc
#

export IMPP_HOME=`pwd`

echo $IMPP_HOME
#Installing FragGeneScan
echo -e "\nInstalling FragGeneScan..."
cd $IMPP_HOME/lib
if [ -d fgs ]; then
	rm -rf fgs
fi
tar xvzf FragGeneScan1.31.tar.gz
cd fgs
make clean
make fgs
[ $? -ne 0 ] && exit $?


## Installing SPAdes
cd $IMPP_HOME/lib
if [ -d spades ]; then
        rm -rf spades
fi
tar xvzf SPAdes-3.15.3-Linux.tar.gz
[ $? -ne 0 ] && exit $?


#Installing BWA
echo -e "\n\nInstalling BWA..."
cd $IMPP_HOME/lib
if [ -d bwa ]; then
        rm -rf bwa
fi
tar xvzf bwa-0.7.17.tar.gz
cd bwa
make
[ $? -ne 0 ] && exit $?


#Installing Plass peptide assembler
echo -e "\n\nInstalling Plass..."
cd $IMPP_HOME/lib
if [ -d plass ]; then
        rm -rf plass
fi
tar xvzf plass-static_sse41.tar.gz

#Installing SGA (binary version)
echo -e "\nInstalling SGA..."
cd $IMPP_HOME/lib
#if [ -f sga]; then
#        rm -f sga
#fi
tar xvzf sga-binary-Linux64.tar.gz

## Installing diamond aligner
echo -e "\n\nInstalling Diamond..."
cd $IMPP_HOME/lib
if [ -d diamond ]; then
        rm -rf diamond
fi
tar xvf diamond.tar.gz

echo -e "\nCompleted installing all third-party tools."
echo -e "\n----------------------------------------------\n"
## IMPP benchmark scripts compile
if [ ! -e $IMPP_HOME/scripts/bin ]; then
        mkdir $IMPP_HOME/scripts/bin
fi
cd $IMPP_HOME/scripts/src
make clean
make all
[ $? -ne 0 ] && exit $?

cd $IMPP_HOME

echo -e "\n-----------------------------------------------\n"
echo -e "\nSuccessfully installed 3rd party tools and benchmarking scripts\n"
echo -e "\n-----------------------------------------------\n"
