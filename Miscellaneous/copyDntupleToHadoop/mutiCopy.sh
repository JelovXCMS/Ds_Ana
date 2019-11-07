#!/bin/sh

INFILELIST=$1
DIR=$2

gfal-mkdir gsiftp://cms-gridftp.rcac.purdue.edu/store/user/chengchi/DsNtuple/${DIR}

while read line
do 
	INPUTFILE=$line
	./copyToHadoop.sh ${INPUTFILE}

done<$INFILELIST
