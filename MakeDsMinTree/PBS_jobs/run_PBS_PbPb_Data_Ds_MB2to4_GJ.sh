#!/bin/sh

INFILELIST=input_list/Ds_pp_Data_MB5.lis
OUTPUTDIR=/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB5
JobIndex=Ds_pp_Data_MB5
FilePerJob=150
isPbPb=3  # 1 is 30-100 trigger
isReal=1

DIR="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree" ## where loop exist
TotalFiles=$(wc -l < "$INFILELIST")
StartFiles=0

echo $TotalFiles , dont use here , not accurate

#cd /home/peng43/work/Project/Ds_PbPb/CMSSW/DsFinder/working_DsFinder/CMSSW_7_5_8_patch3/src
#eval `scramv1 runtime -sh`
CMSenv
cd ${DIR}/PBS_jobs

#rm ../MakeDsMinTree.exe
#g++ ../MakeDsMinTree.C $(root-config --cflags --libs) -Wall -O2 -o ../MakeDsMinTree.exe


cd $DIR/PBS_jobs
#./submit_PBS_job_mutltifiles.sh INFILELIST isPbPb FilePerJob OUTPUTDIR JobIndex isReal
#../submit_PBS_job_mutltifiles.sh $INFILELIST $isPbPb $FilePerJob $OUTPUTDIR $JobIndex $isReal

../submit_PBS_job_mutltifiles.sh input_list/Ds_PbPb_Data_HIMB2_GJ.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/GJ/Ds_PbPb_Data_HIMB2_GJ Ds_PbPb_Data_HIMB2_GJ  $isReal
../submit_PBS_job_mutltifiles.sh input_list/Ds_PbPb_Data_HIMB3_GJ.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/GJ/Ds_PbPb_Data_HIMB3_GJ Ds_PbPb_Data_HIMB3_GJ  $isReal
../submit_PBS_job_mutltifiles.sh input_list/Ds_PbPb_Data_HIMB4_GJ.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/GJ/Ds_PbPb_Data_HIMB4_GJ Ds_PbPb_Data_HIMB4_GJ  $isReal




## ./exec_condorfilelist.sh		inputFilelist		isPbPb		startfile		endfile		outputfolder 	isReal	 
##./exec_condorfilelist.sh $INFILELIST $isPbPb $StartFiles $TotalFiles $OUTPUTDIR $isReal

cd ${DIR}/PBS_jobs 
