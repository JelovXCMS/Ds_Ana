#!/bin/bash


#cd /home/peng43/root_from_CMSSW/root602_13/CMSSW_7_5_8_patch3/src
#eval `scramv1 runtime -sh`
#cd -


cd /home/peng43/root_from_CMSSW/root610_09/CMSSW_10_1_0_pre1/src
eval `scramv1 runtime -sh`
cd -

rm Fit_sideband.exe
# g++ Fit_sideband.C $(root-config --cflags --libs) -L $ROOTSYS/lib -lRooFit -lHtml -lMinuit -lRooFitCore -Wall -O2 -o Fit_sideband.C.exe
# g++ $(root-config --cflags --libs) -L $ROOTSYS/lib -lRooFit -lHtml -lMinuit -lRooFitCore Fit_sideband.C #work
g++ Fit_sideband.C $(root-config --cflags --libs) -L $ROOTSYS/lib -lRooFit -lHtml -lMinuit -lRooFitCore -Wall -O2 -o Fit_sideband.exe
# g++ $(root-config --cflags --libs) -L $ROOTSYS/lib -lRooFit -lHtml -lMinuit -lRooFitCore c.cpp
# g++ Fit_sideband.C $(root-config --cflags --libs)-L $ROOTSYS/lib -lRooFit -lHtml -lMinuit  -lRooFitCore -o Fit_sideband.exe



