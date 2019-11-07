#!/bin/sh

#int MC_MergeGenSamplePtWeight(int isPbPb=1, int PNPrompt =0, int DsChannel=0, TString MCFileList="./MC_List/DsMinTree_PbPb_MC_Prompt_phikkpi.lis", TString inMergeFile="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/PbPb_MC/DsMinTree_PbPb_MC_Prompt_phikkpi.root")


# root -b -l -q 'MC_MergeGenSamplePtWeight.C++ (0,0,0,"./MC_List/DsMinTree_pp_offMC_Prompt_phikkpi.lis","/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/output/DsMinTree_pp_offMC_Ds_Prompt_phikkpi_Merge.root")'      > log/pp_Prompt_phi.log

mkdir -p /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_offMC_nocut/

# pp
root -b -l -q 'MC_MergeGenSamplePtWeight.C++ (0,0,0,"./MC_List/DsMinTree_pp_offMC_Prompt_Phi_nocut.lis","/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_offMC_nocut/DsMinTree_pp_offMC_Prompt_phi.root")'      > log/pp_Prompt_phi_nocut.log 
##root -b -l -q 'MC_MergeGenSamplePtWeight.C++ (0,0,1,"./MC_List/DsMinTree_pp_offMC_Prompt_f0.lis","/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_offMC/DsMinTree_pp_offMC_Prompt_f0.root")'      > log/pp_Prompt_f0.log 
root -b -l -q 'MC_MergeGenSamplePtWeight.C++ (0,1,0,"./MC_List/DsMinTree_pp_offMC_NonPrompt_Phi_nocut.lis","/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_offMC_nocut/DsMinTree_pp_offMC_NonPrompt_phi.root")'      > log/pp_NonPrompt_phi_nocut.log 
##root -b -l -q 'MC_MergeGenSamplePtWeight.C++ (0,1,1,"./MC_List/DsMinTree_pp_offMC_NonPrompt_f0.lis","/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_offMC/DsMinTree_pp_offMC_NonPrompt_f0.root")'      > log/pp_NonPrompt_f0.log 

wait

#echo "pp done"

## PbPb
root -b -l -q 'MC_MergeGenSamplePtWeight.C+ (1,0,0,"./MC_List/DsMinTree_PbPb_offMC_Prompt_Phi_nocut.lis","/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_offMC_nocut/DsMinTree_PbPb_offMC_Prompt_phi.root")'      > log/PbPb_Prompt_phi_nocut.log &
#root -b -l -q 'MC_MergeGenSamplePtWeight.C+ (1,0,1,"./MC_List/DsMinTree_PbPb_offMC_Prompt_f0.lis","/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_offMC/DsMinTree_PbPb_offMC_Prompt_f0.root")'      > log/PbPb_Prompt_f0.log &
#
root -b -l -q 'MC_MergeGenSamplePtWeight.C+ (1,1,0,"./MC_List/DsMinTree_PbPb_offMC_NonPrompt_Phi_nocut.lis","/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_offMC_nocut/DsMinTree_PbPb_offMC_NonPrompt_phi.root")'      > log/PbPb_NonPrompt_phi_nocut.log &
#root -b -l -q 'MC_MergeGenSamplePtWeight.C+ (1,1,1,"./MC_List/DsMinTree_PbPb_offMC_NonPrompt_f0.lis","/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_offMC/DsMinTree_PbPb_offMC_NonPrompt_f0.root")'      > log/PbPb_NonPrompt_f0.log &
#
wait


#root -b -l -q 'MC_MergeGenSamplePtWeight.C++ (0,0,0,"./MC_List/DsMinTree_pp_MC_Prompt_phikkpi.lis")'      > log/pp_Prompt_phi.log
#root -b -l -q 'MC_MergeGenSamplePtWeight.C++ (0,0,1,"./MC_List/DsMinTree_pp_MC_Prompt_f0980kkpi.lis")'    > log/pp_Prompt_f0980.log
#root -b -l -q 'MC_MergeGenSamplePtWeight.C++ (0,1,0,"./MC_List/DsMinTree_pp_MC_NonPrompt_phikkpi.lis")'   > log/pp_NonPrompt_phi.log
#root -b -l -q 'MC_MergeGenSamplePtWeight.C++ (0,1,1,"./MC_List/DsMinTree_pp_MC_NonPrompt_f0980kkpi.lis")' > log/pp_NonPrompt_f0980.log
#
#echo "pp done"
#
#root -b -l -q 'MC_MergeGenSamplePtWeight.C++ (1,0,0,"./MC_List/DsMinTree_PbPb_MC_Prompt_phikkpi.lis")'      > log/PbPb_Prompt_phi.log
#
#echo "PbPb Prompt phi done"
#
#root -b -l -q 'MC_MergeGenSamplePtWeight.C++ (1,0,1,"./MC_List/DsMinTree_PbPb_MC_Prompt_f0980kkpi.lis")'    > log/PbPb_Prompt_f0980.log
#
#echo "PbPb Prompt f0 done"
#
#root -b -l -q 'MC_MergeGenSamplePtWeight.C++ (1,1,0,"./MC_List/DsMinTree_PbPb_MC_NonPrompt_phikkpi.lis")'   > log/PbPb_NonPrompt_phi.log
#
#echo "PbPb NonPrompt phi done"
#
#root -b -l -q 'MC_MergeGenSamplePtWeight.C++ (1,1,1,"./MC_List/DsMinTree_PbPb_MC_NonPrompt_f0980kkpi.lis")' > log/PbPb_NonPrompt_f0980.log
#
#echo "PbPb NonPrompt f0 done"
