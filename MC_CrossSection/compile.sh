#!/bin/sh


cd /home/peng43/root_from_CMSSW/root602_13/CMSSW_7_5_8_patch3/src
eval `scramv1 runtime -sh`
cd -

rm CSProject.exe
g++ CSProject.C $(root-config --cflags --libs) -Wall -O2 -o CSProject.exe


# rm MakeDsMinTreeV2.exe
# g++ MakeDsMinTreeV2.C $(root-config --cflags --libs) -Wall -O2 -o MakeDsMinTreeV2.exe

