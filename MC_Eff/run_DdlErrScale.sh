#!/bin/sh
#int MC_EffCalLoop(int isPbPb=0, int PNPrompt = 0, int DsChannel=0  ){


root -b -q 'MC_EffCalLoop_DdlErrScale.C++ (0,0,0)' >log/pp_Prompt_phi.log  
root -b -q 'MC_EffCalLoop_DdlErrScale.C++ (0,1,0)' >log/pp_NonPrompt_phi.log    
# root -b -q 'MC_EffCalLoop_DdlErrScale.C++ (0,0,1)' >log/pp_Prompt_f0980.log  
# root -b -q 'MC_EffCalLoop_DdlErrScale.C++ (0,1,1)' >log/pp_NonPrompt_f0980.log   

echo "pp done"

root -b -q 'MC_EffCalLoop_DdlErrScale.C++ (3,0,0)' >log/PbPb_Prompt_phi.log  
root -b -q 'MC_EffCalLoop_DdlErrScale.C++ (3,1,0)' >log/PbPb_NonPrompt_phi.log   
#echo "PbPb phi done"
#root -b -q 'MC_EffCalLoop_DdlErrScale.C++ (3,0,1)' >log/PbPb_Prompt_f0980.log  
#root -b -q 'MC_EffCalLoop_DdlErrScale.C++ (3,1,1)' >log/PbPb_NonPrompt_f0980.log  
echo "PbPb f0 done"


