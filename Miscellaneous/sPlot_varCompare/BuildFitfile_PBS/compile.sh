#!/bin/sh

cd /home/peng43/root_from_CMSSW/root602_13/CMSSW_7_5_8_patch3/src
#cd /home/peng43/root_from_CMSSW/root610_09/CMSSW_10_1_0_pre1/src
eval `scramv1 runtime -sh`
cd -

rm BuildFitfile.exe
g++ BuildFitfile.C $(root-config --cflags --libs) -Wall -O2 -o BuildFitfile.exe

