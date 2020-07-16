#bin/bash
#
# Script to install third-party softwares
#------------------------
# Pre-requisites needed
#------------------------
# cmake :
# zlib :
# gcc version:
#
#
#

export IMPP_HOME=`pwd`

echo $IMPP_HOME
#Installing FragGeneScan
echo -e "\nInstalling FragGeneScan..."
cd $IMPP_HOME/lib
if [ -d FragGeneScan1.31 ]; then
	rm -rf FragGeneScan1.31
fi
tar xvzf FragGeneScan1.31.tar.gz
cd FragGeneScan1.31
make clean
make fgs
[ $? -ne 0 ] && exit $?


#Installing SPAdes
echo -e "\n\nInstalling SPAdes..."
cd $IMPP_HOME/lib
if [ -d SPAdes-3.9.0-Linux ]; then
        rm -rf SPAdes-3.9.0-Linux
fi
tar xvzf SPAdes-3.9.0-Linux.tar.gz
[ $? -ne 0 ] && exit $?


#Installing BWA
echo -e "\n\nInstalling BWA..."
cd $IMPP_HOME/lib
if [ -d bwa-0.7.17 ]; then
        rm -rf bwa-0.7.17
fi
tar xvf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
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


echo -e "\nCompleted installing all third-party tools."
echo -e "\n----------------------------------------------\n"
## IMPP installation
echo -e "\nInstalling iMPP ..."
if [ ! -e $IMPP_HOME/bin ]; then
        mkdir $IMPP_HOME/bin
fi
cd $IMPP_HOME/src
make clean
make all
[ $? -ne 0 ] && exit $?

cd $IMPP_HOME

echo -e "\n-----------------------------------------------\n"
echo -e "\nSuccessfully installed iMPP.\n"
echo -e "\n-----------------------------------------------\n"
