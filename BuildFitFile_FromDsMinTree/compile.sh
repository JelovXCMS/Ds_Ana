#!/bin/sh

cd /home/peng43/root_from_CMSSW/root602_13/CMSSW_7_5_8_patch3/src
#cd /home/peng43/root_from_CMSSW/root610_09/CMSSW_10_1_0_pre1/src
eval `scramv1 runtime -sh`
cd -

rm BuildFitFile_FromDsMinTreeLoop.exe
g++ BuildFitFile_FromDsMinTreeLoop.C $(root-config --cflags --libs) -Wall -O2 -o BuildFitFile_FromDsMinTreeLoop.exe

rm BuildFitFile_FromDsMinTreeLoop_ForCut.exe
g++ BuildFitFile_FromDsMinTreeLoop_ForCut.C $(root-config --cflags --libs) -Wall -O2 -o BuildFitFile_FromDsMinTreeLoop_ForCut.exe


rm BuildFitFile_FromDsMinTreeLoop_moreScan.exe
g++ BuildFitFile_FromDsMinTreeLoop_moreScan.C $(root-config --cflags --libs) -Wall -O2 -o BuildFitFile_FromDsMinTreeLoop_moreScan.exe


