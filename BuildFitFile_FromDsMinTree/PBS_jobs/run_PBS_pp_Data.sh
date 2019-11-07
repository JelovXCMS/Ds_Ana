#!/bin/sh

INFILELIST=input_list/pp_Data_DsMinTree.lis
OUTPUTDIR=/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BuildFitFile_FromDsMinTree/output
JobIndex=BuildFit_pp_GJ
FilePerJob=1
isPbPb=0  # 1 is 30-100 trigger
isReal=1

DIR="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BuildFitFile_FromDsMinTree" ## where loop exist
TotalFiles=$(wc -l < "$INFILELIST")
StartFiles=0

echo $TotalFiles , dont use here , not accurate

#cd /home/peng43/work/Project/Ds_PbPb/CMSSW/DsFinder/working_DsFinder/CMSSW_7_5_8_patch3/src
#eval `scramv1 runtime -sh`
#CMSenvNew
cd ${DIR}/PBS_jobs

#rm ../MakeDsMinTree.exe
#g++ ../MakeDsMinTree.C $(root-config --cflags --libs) -Wall -O2 -o ../MakeDsMinTree.exe


cd $DIR/PBS_jobs
#./submit_PBS_job_mutltifiles.sh INFILELIST isPbPb FilePerJob OUTPUTDIR JobIndex isReal
#../submit_PBS_job_mutltifiles.sh $INFILELIST $isPbPb $FilePerJob $OUTPUTDIR $JobIndex $isReal

../submit_PBS_job_mutltifiles.sh $INFILELIST  $isPbPb $FilePerJob $OUTPUTDIR $JobIndex  $isReal

## ./exec_condorfilelist.sh		inputFilelist		isPbPb		startfile		endfile		outputfolder 	isReal	 
##./exec_condorfilelist.sh $INFILELIST $isPbPb $StartFiles $TotalFiles $OUTPUTDIR $isReal

cd ${DIR}/PBS_jobs 