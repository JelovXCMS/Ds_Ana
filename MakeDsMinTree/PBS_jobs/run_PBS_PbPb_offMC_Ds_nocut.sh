#!/bin/sh

INFILELIST=input_list/Ds_pp_Data_MB5_nocut.lis
OUTPUTDIR=/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB5
JobIndex=Ds_pp_Data_MB5
FilePerJob=100
isPbPb=1
isReal=0

DIR="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree" ## where loop exist
TotalFiles=$(wc -l < "$INFILELIST")
StartFiles=0


ptbins=(0 2 4 6 10 19)

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

for i in {0..5}
do
	../submit_PBS_job_mutltifiles.sh input_list/Ds_PbPb_MC_Prompt_Phi_pt${ptbins[i]}_nocut.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_offMC_nocut/Ds_PbPb_MC_Prompt_Phi_pt${ptbins[i]} Ds_PbPb_MC_Prompt_Phi_pt${ptbins[i]}  $isReal

	../submit_PBS_job_mutltifiles.sh input_list/Ds_PbPb_MC_NonPrompt_Phi_pt${ptbins[i]}_nocut.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_offMC_nocut/Ds_PbPb_MC_NonPrompt_Phi_pt${ptbins[i]} Ds_PbPb_MC_NonPrompt_Phi_pt${ptbins[i]}  $isReal

#	../submit_PBS_job_mutltifiles.sh input_list/Ds_PbPb_MC_Prompt_f0_pt${ptbins[i]}_nocut.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_offMC_nocut/Ds_PbPb_MC_Prompt_f0_pt${ptbins[i]} Ds_PbPb_MC_Prompt_f0_pt${ptbins[i]}  $isReal
 
#	../submit_PBS_job_mutltifiles.sh input_list/Ds_PbPb_MC_NonPrompt_f0_pt${ptbins[i]}_nocut.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_offMC_nocut/Ds_PbPb_MC_NonPrompt_f0_pt${ptbins[i]} Ds_PbPb_MC_NonPrompt_f0_pt${ptbins[i]}  $isReal

done

#../submit_PBS_job_mutltifiles.sh input_list/Ds_PbPb_Data_HIMB5_nocut.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/cent30100/Ds_PbPb_Data_HIMB5 Ds_PbPb_Data_HIMB5  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_PbPb_Data_HIMB6_nocut.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/cent30100/Ds_PbPb_Data_HIMB6 Ds_PbPb_Data_HIMB6  $isReal
#../submit_PBS_job_mutltifiles.sh input_list/Ds_PbPb_Data_HIMB7_nocut.lis  $isPbPb $FilePerJob /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/cent30100/Ds_PbPb_Data_HIMB7 Ds_PbPb_Data_HIMB7  $isReal


## ./exec_condorfilelist.sh		inputFilelist		isPbPb		startfile		endfile		outputfolder 	isReal	 
##./exec_condorfilelist.sh $INFILELIST $isPbPb $StartFiles $TotalFiles $OUTPUTDIR $isReal

cd ${DIR}/PBS_jobs 
