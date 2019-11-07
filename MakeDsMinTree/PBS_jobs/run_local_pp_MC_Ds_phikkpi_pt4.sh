#!/bin/sh

INFILELIST=input_list/Ds_phikkpi_pt4.lis
OUTPUTDIR=/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_MC/Ds_phikkpi_pt4
TotalFiles=$(wc -l < "$INFILELIST")
StartFiles=0
isPbPb=0
isReal=0
DIR="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree" ## where loop exist

echo $TotalFiles

#eval `scramv1 runtime -sh`
CMSenv
#rm ../MakeDsMinTree.exe
#g++ ../MakeDsMinTree.C $(root-config --cflags --libs) -Wall -O2 -o ../MakeDsMinTree.exe


cd $DIR

## ./exec_condorfilelist.sh		inputFilelist		isPbPb		startfile		endfile		outputfolder 	isReal	 
./exec_condorfilelist.sh $INFILELIST $isPbPb $StartFiles $TotalFiles $OUTPUTDIR $isReal

cd PBS_jobs 
