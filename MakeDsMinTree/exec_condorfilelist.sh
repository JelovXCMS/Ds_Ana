#!/bin/sh

echo `hostname`
echo `date`
export SCRAM_ARCH=slc6_amd64_gcc491
# source /apps/osg/cmssoft/cmsset_default.sh
source /cvmfs/cms.cern.ch/cmsset_default.sh 
export X509_USER_PROXY=/home/peng43/.myproxy
DIR="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree"   ## this is where loop.exe exist
cd /home/peng43/root_from_CMSSW/root602_13/CMSSW_7_5_8_patch3/src
eval `scramv1 runtime -sh`
cd $DIR

cd /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree

INFILELIST=$1
isPbPb=$2
STARTFILE=$3
ENDFILE=$4
OUTPUTDIR=$5
isReal=$6

INFILELIST=${DIR}/PBS_jobs/${INFILELIST}

echo "INFILELIST: $INFILELIST"
echo "STARTFILE: $STARTFILE, ENDFILE: $ENDFILE"

if [ ! -d $OUTPUTDIR ]; then
	mkdir -p $OUTPUTDIR
fi

fileCounter=-1

while read line
do
	fileCounter=$((fileCounter+1))
	if [ $fileCounter -lt $STARTFILE ] || [ $fileCounter -ge $ENDFILE ]; then
		continue
	fi
	echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
	echo "<<<<<<<<<<<<< NEW INPUTFILE >>>>>>>>>>>>>>>>>"
	echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
	echo "fileCounter: $fileCounter, inputfile: $line"

	Inputfile=$line
	Outputfile="$OUTPUTDIR/DsMinTree_$(basename ${Inputfile})_${fileCounter}.root"
	echo "Outputfile: $Outputfile"

	./MakeDsMinTree.exe $Inputfile $Outputfile $isReal $isPbPb

	#./something.exe for other purpose
	
	echo "<<<<<<<<<<<<< DONE!!!!!!! <<<<<<<<<<<<<<<<<"
done<$INFILELIST
echo `date`
