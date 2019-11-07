#!/bin/sh

INFILELIST=input_list/Ds_pp_Data_MB5.lis
OUTPUTDIR=/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB5
JobIndex=Ds_pp_Data_MB5
FilePerJob=10
isPbPb=0
isReal=1

DIR="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_CrossSection" ## where loop exist
TotalFiles=$(wc -l < "$INFILELIST")
StartFiles=0

echo $TotalFiles

cd /home/peng43/work/Project/Ds_PbPb/CMSSW/DsFinder/working_DsFinder/CMSSW_7_5_8_patch3/src
eval `scramv1 runtime -sh`
#CMSenv
cd ${DIR}/PBS_jobs

#rm ../MakeDsMinTree.exe
#g++ ../MakeDsMinTree.C $(root-config --cflags --libs) -Wall -O2 -o ../MakeDsMinTree.exe


cd $DIR/PBS_jobs
#./submit_PBS_job_mutltifiles.sh INFILELIST isPbPb FilePerJob OUTPUTDIR JobIndex isReal
#../submit_PBS_job_mutltifiles.sh $INFILELIST $isPbPb $FilePerJob $OUTPUTDIR $JobIndex $isReal

../submit_PBS_job_mutltifiles.sh input_list/CMStune_CharmAll_pthat20.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/CS_CharmAll_ph20 CS_wholeTable_CharmAllph20  $isReal

## ./exec_condorfilelist.sh		inputFilelist		isPbPb		startfile		endfile		outputfolder 	isReal	 
##./exec_condorfilelist.sh $INFILELIST $isPbPb $StartFiles $TotalFiles $OUTPUTDIR $isReal

cd ${DIR}/PBS_jobs 
