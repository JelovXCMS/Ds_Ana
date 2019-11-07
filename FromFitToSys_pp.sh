#!/bin/sh


cd ./SignalFit
root -b -q ResultSum.C++(0)
cd ..
cd ./CrossSection_dNdpt/
root -b -q CSdNdpt.C++(0)
cd ..
cd ./Systematics/
./run_pp.sh
#root -b -q CutScanSys.C++(0) 
cd ..
