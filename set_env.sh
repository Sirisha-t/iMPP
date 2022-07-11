#!/bin/bash

ulimit -s unlimited

export IMPP_HOME=`pwd`

if [[ $IMPP_HOME == "" ]]; then
        echo "IMPP_HOME is not properly configured"
        exit
fi

export IMPP_BIN=$IMPP_HOME/bin
export IMPP_SCRIPT=$IMPP_HOME/wrapper
export FGS=$IMPP_HOME/lib/FragGeneScan1.31
export PLASS=$IMPP_HOME/lib/plass/bin/plass
export BWA=$IMPP_HOME/bwa-0.7.17/bwa
export SGA=$IMPP_HOME/lib/sga
export SPADES=$IMPP_HOME/SPAdes-3.9.0-Linux/bin/metaspades.py
export PATH=$IMPP_HOME:$IMPP_BIN:$IMPP_SCRIPT:$FGS:$PLASS:$BWA:$SGA:$SPADES:$PATH

