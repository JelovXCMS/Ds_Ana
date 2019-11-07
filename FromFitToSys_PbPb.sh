#!/bin/sh


cd ./SignalFit
root -b -q ResultSum.C++(3)
cd ..
cd ./CrossSection_dNdpt/
root -b -q CSdNdpt.C++(3)
cd ..
cd ./Systematics/
#root -b -q CutScanSys.C++(3) 
./run_PbPb.sh
cd ..
