#!/bin/sh

# int BuildFitfile(Int_t isPbPb=0, Int_t isReal=0,int PNPrompt=0,int DsChannel=0, TString s_fin="/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/DsMinTree_pp_Data_MBAll.root",TString foutTag="", double Dpt_Low=0, double Dpt_High=0)

mkdir -p log

root -b -q 'BuildFitfile.C++(0,1,0,0,"/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/DsMinTree_pp_Data_MBAll.root","")' > log/BF_pp_Data.log 

echo "pp Data done"

# root -b -q 'BuildFitfile.C++(0,0,0,0)' > log/BF_pp_MC_Prompt_phi.log
# root -b -q 'BuildFitfile.C++(0,0,0,1)' > log/BF_pp_MC_Prompt_f0.log

echo "pp MC prompt done"

# root -b -q 'BuildFitfile.C++(0,0,1,0)' > log/BF_pp_MC_NonPrompt_phi.log
# root -b -q 'BuildFitfile.C++(0,0,1,1)' > log/BF_pp_MC_NonPrompt_f0.log

echo "pp MC nonprompt done"
