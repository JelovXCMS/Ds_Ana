#!/bin/sh


cd /home/peng43/root_from_CMSSW/root602_13/CMSSW_7_5_8_patch3/src
eval `scramv1 runtime -sh`
cd -

rm MakeDsMinTree.exe
g++ MakeDsMinTree.C $(root-config --cflags --libs) -Wall -O2 -o MakeDsMinTree.exe


rm MakeDsMinTreeV2.exe
g++ MakeDsMinTreeV2.C $(root-config --cflags --libs) -Wall -O2 -o MakeDsMinTreeV2.exe

