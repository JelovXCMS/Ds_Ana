#!/bin/sh

INFILELIST=input_list/Ds_pp_Data_MB5.lis
OUTPUTDIR=/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB5
JobIndex=Ds_pp_Data_MB5
FilePerJob=100
isPbPb=0
isReal=1

DIR="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree" ## where loop exist
#TotalFiles=$(wc -l < "$INFILELIST")
StartFiles=0

#echo $TotalFiles

cd /home/peng43/work/Project/Ds_PbPb/CMSSW/DsFinder/working_DsFinder/CMSSW_7_5_8_patch3/src
eval `scramv1 runtime -sh`
#CMSenv
cd ${DIR}/PBS_jobs

#rm ../MakeDsMinTree.exe
#g++ ../MakeDsMinTree.C $(root-config --cflags --libs) -Wall -O2 -o ../MakeDsMinTree.exe


cd $DIR/PBS_jobs
#./submit_PBS_job_mutltifiles.sh INFILELIST isPbPb FilePerJob OUTPUTDIR JobIndex isReal
#../submit_PBS_job_mutltifiles.sh $INFILELIST $isPbPb $FilePerJob $OUTPUTDIR $JobIndex $isReal

for i in {1..20}
do
../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB${i}.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB${i} Ds_pp_Data_MB${i}  $isReal
done

#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB2.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB2 Ds_pp_Data_MB2  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB3.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB3 Ds_pp_Data_MB3  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB4.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB4 Ds_pp_Data_MB4  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB5.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB5 Ds_pp_Data_MB5  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB6.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB6 Ds_pp_Data_MB6  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB7.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB7 Ds_pp_Data_MB7  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB8.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB8 Ds_pp_Data_MB8  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB9.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB9 Ds_pp_Data_MB9  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB10.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB10 Ds_pp_Data_MB10  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB11.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB11 Ds_pp_Data_MB11  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB12.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB12 Ds_pp_Data_MB12  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB13.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB13 Ds_pp_Data_MB13  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB14.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB14 Ds_pp_Data_MB14  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB15.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB15 Ds_pp_Data_MB15  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB16.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB16 Ds_pp_Data_MB16  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB17.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB17 Ds_pp_Data_MB17  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB18.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB18 Ds_pp_Data_MB18  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB19.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB19 Ds_pp_Data_MB19  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_pp_Data_MB20.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB20 Ds_pp_Data_MB20  $isReal


## ./exec_condorfilelist.sh		inputFilelist		isPbPb		startfile		endfile		outputfolder 	isReal	 
##./exec_condorfilelist.sh $INFILELIST $isPbPb $StartFiles $TotalFiles $OUTPUTDIR $isReal

cd ${DIR}/PBS_jobs 
