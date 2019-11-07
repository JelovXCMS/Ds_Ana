#!/bin/bash


isPbPb=0
DptLow=6
DptHigh=40
# DvtxCut=0.05
# DdlsCut=2.5

DVtxPCut=0.25
DdlsCut=4.5


doProject=1
doMkkpiFit=1
doMkkFit=1

if [$doProject -eq 1]; then
root -b -q ProjectToHis.C++(isPbPb,DptLow,DptHigh,DvtxCut,DdlsCut)

fi

# if [$doMkkpiFit -eq 1]; then


# fi

# if [$doMkkFit -eq 1]; then


# fi

