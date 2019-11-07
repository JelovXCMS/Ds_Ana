#!/bin/sh

#int plot_weight(int isPbPb=1, int PNPrompt =0, int DsChannel=0, TString MCFileList="./MC_List/DsMinTree_PbPb_MC_Prompt_phikkpi.lis", TString inMergeFile="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/PbPb_MC/DsMinTree_PbPb_MC_Prompt_phikkpi.root")


root -b -l -q 'plot_weight.C++ (0,0,0,"./MC_List/DsMinTree_pp_offMC_Prompt_phikkpi.lis")'      
#root -b -l -q 'plot_weight.C++ (0,0,1,"./MC_List/DsMinTree_pp_offMC_Prompt_f0980kkpi.lis")'    
root -b -l -q 'plot_weight.C++ (0,1,0,"./MC_List/DsMinTree_pp_offMC_NonPrompt_phikkpi.lis")'   
#root -b -l -q 'plot_weight.C++ (0,1,1,"./MC_List/DsMinTree_pp_offMC_NonPrompt_f0980kkpi.lis")' 

echo "pp done"

root -b -l -q 'plot_weight.C++ (1,0,0,"./MC_List/DsMinTree_PbPb_offMC_Prompt_phikkpi.lis")'   
#root -b -l -q 'plot_weight.C++ (1,0,1,"./MC_List/DsMinTree_PbPb_offMC_Prompt_f0980kkpi.lis")'  
root -b -l -q 'plot_weight.C++ (1,1,0,"./MC_List/DsMinTree_PbPb_offMC_NonPrompt_phikkpi.lis")'  
#root -b -l -q 'plot_weight.C++ (1,1,1,"./MC_List/DsMinTree_PbPb_offMC_NonPrompt_f0980kkpi.lis")' 

echo "PbPb NonPrompt f0 done"
