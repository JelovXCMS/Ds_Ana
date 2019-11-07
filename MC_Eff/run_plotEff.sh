#!/bin/sh
#int MC_EffCalLoop(int isPbPb=0, int PNPrompt = 0, int DsChannel=0  ){


root -b -q 'Plot_Eff.C++ (0,0)' 
# root -b -q 'Plot_Eff.C++ (0,1)' // this is for f0 
root -b -q 'Plot_Eff.C++ (3,0)' 
# root -b -q 'Plot_Eff.C++ (3,1)' 
echo "PbPb f0 done"


