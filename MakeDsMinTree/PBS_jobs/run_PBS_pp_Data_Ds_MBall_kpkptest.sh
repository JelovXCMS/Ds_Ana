#!/bin/sh

#INFILELIST=input_list/Ds_pp_Data_MB5.lis
#OUTPUTDIR=/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB5
#JobIndex=Ds_pp_Data_MB5
FilePerJob=40
isPbPb=0
isReal=1

DIR="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree" ## where loop exist
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

../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB1_kpkptest.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/Ds_pp_Data_MB1_kpkptest Ds_pp_Data_MB1  $isReal

#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB2_kpkptest.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/Ds_pp_Data_MB2_kpkptest Ds_pp_Data_MB2  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB3_kpkptest.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/Ds_pp_Data_MB3_kpkptest Ds_pp_Data_MB3  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB4_kpkptest.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/Ds_pp_Data_MB4_kpkptest Ds_pp_Data_MB4  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB5_kpkptest.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/Ds_pp_Data_MB5_kpkptest Ds_pp_Data_MB5  $isReal

../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB6_kpkptest.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/Ds_pp_Data_MB6_kpkptest Ds_pp_Data_MB6  $isReal
../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB7_kpkptest.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/Ds_pp_Data_MB7_kpkptest Ds_pp_Data_MB7  $isReal
../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB8_kpkptest.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/Ds_pp_Data_MB8_kpkptest Ds_pp_Data_MB8  $isReal
../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB9_kpkptest.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/Ds_pp_Data_MB9_kpkptest Ds_pp_Data_MB9  $isReal
../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB10_kpkptest.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/Ds_pp_Data_MB10_kpkptest Ds_pp_Data_MB10  $isReal
../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB11_kpkptest.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/Ds_pp_Data_MB11_kpkptest Ds_pp_Data_MB11  $isReal
../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB12_kpkptest.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/Ds_pp_Data_MB12_kpkptest Ds_pp_Data_MB12  $isReal
../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB13_kpkptest.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/Ds_pp_Data_MB13_kpkptest Ds_pp_Data_MB13  $isReal
../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB14_kpkptest.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/Ds_pp_Data_MB14_kpkptest Ds_pp_Data_MB14  $isReal
../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB15_kpkptest.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/Ds_pp_Data_MB15_kpkptest Ds_pp_Data_MB15  $isReal
../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB16_kpkptest.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/Ds_pp_Data_MB16_kpkptest Ds_pp_Data_MB16  $isReal
../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB17_kpkptest.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/Ds_pp_Data_MB17_kpkptest Ds_pp_Data_MB17  $isReal
../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB18_kpkptest.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/Ds_pp_Data_MB18_kpkptest Ds_pp_Data_MB18  $isReal
../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB19_kpkptest.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/Ds_pp_Data_MB19_kpkptest Ds_pp_Data_MB19  $isReal
../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB20_kpkptest.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/Ds_pp_Data_MB20_kpkptest Ds_pp_Data_MB20  $isReal

## ./exec_condorfilelist.sh		inputFilelist		isPbPb		startfile		endfile		outputfolder 	isReal	 
##./exec_condorfilelist.sh $INFILELIST $isPbPb $StartFiles $TotalFiles $OUTPUTDIR $isReal

cd ${DIR}/PBS_jobs 
