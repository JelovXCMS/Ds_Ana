#!/bin/sh

INFILELIST=input_list/PbPb_Data_DsMinTree.lis
OUTPUTDIR=/scratch/halstead/p/peng43/Ds_phikkpi/BuildFitFile_ForVarCompare/PbPb3_Data/
JobIndex=BuildFitVarCompare_PbPb3
FilePerJob=400
isPbPb=3  # 1 is 30-100 trigger
isReal=1

DIR="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/Miscellaneous/sPlot_varCompare/BuildFitfile_PBS" ## where loop exist
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

#../submit_PBS_job_mutltifiles.sh $INFILELIST  $isPbPb $FilePerJob $OUTPUTDIR $JobIndex  $isReal

../submit_PBS_job_mutltifiles.sh input_list/Ds_PbPb_Data_HIMB2_GJ.lis       $isPbPb $FilePerJob ${OUTPUTDIR}/HIMB2_GJ/ HIMB2GJ${JobIndex} $isReal
../submit_PBS_job_mutltifiles.sh input_list/Ds_PbPb_Data_HIMB3_GJ.lis       $isPbPb $FilePerJob ${OUTPUTDIR}/HIMB3_GJ/ HIMB3GJ${JobIndex} $isReal
../submit_PBS_job_mutltifiles.sh input_list/Ds_PbPb_Data_HIMB4_GJ.lis       $isPbPb $FilePerJob ${OUTPUTDIR}/HIMB4_GJ/ HIMB4GJ${JobIndex} $isReal

for i in {1..11}
	do
	../submit_PBS_job_mutltifiles.sh input_list/Ds_PbPb_Data_HIMB${i}_TrkOnly.lis  $isPbPb $FilePerJob $OUTPUTDIR/HIMB${i}_TrkOnly/ HIMBi${i}TrkOnly$JobIndex  $isReal
done

## ./exec_condorfilelist.sh		inputFilelist		isPbPb		startfile		endfile		outputfolder 	isReal	 
##./exec_condorfilelist.sh $INFILELIST $isPbPb $StartFiles $TotalFiles $OUTPUTDIR $isReal

cd ${DIR}/PBS_jobs 
