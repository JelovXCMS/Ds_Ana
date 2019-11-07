#!/bin/sh

# int BuildFitfile(Int_t isPbPb=0, Int_t isReal=0,int PNPrompt=0,int DsChannel=0, TString s_fin="/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/DsMinTree_pp_Data_MBAll.root",TString foutTag="", double Dpt_Low=0, double Dpt_High=0)

mkdir -p log

#root -b -q 'BuildFitfile.C++(3,1,0,0,"/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/Ds_PbPb_Data_HIMBAll.root","")' > log/BF_PbPb_Data.log & 

sleep 30

echo "PbPb Data done"

root -b -q 'BuildFitfile.C++(3,0,0,0)' > log/BF_PbPb_MC_Prompt_phi.log 
# root -b -q 'BuildFitfile.C+(3,0,0,1)' > log/BF_PbPb_MC_Prompt_f0.log &

wait
echo "PbPb MC prompt done"

root -b -q 'BuildFitfile.C++(3,0,1,0)' > log/BF_PbPb_MC_NonPrompt_phi.log 
# root -b -q 'BuildFitfile.C+(3,0,1,1)' > log/BF_PbPb_MC_NonPrompt_f0.log &

echo "PbPb MC nonprompt done"


wait
