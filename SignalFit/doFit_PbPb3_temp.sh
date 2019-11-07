#!/bin/bash

CMSenvNew

PTBIN=(4 5 6 8 10 20 40)
RAA=(0.60 0.43 0.36 0.27 0.26 0.27 0.32 0.30 0.38)

#COLSYST=(0)  # isPbPb , 0 :pp 1:30-100% , 2:30-100% 3:0-100%

nPT=$((${#PTBIN[@]}-1))
nCOL=${#COLSYST[@]}

i=2
while ((i<$nPT))
#while ((i<4))
do
echo "Processing pT bins: ${PTBIN[i]} - ${PTBIN[i+1]} GeV/c"
#root -l -b -q 'SignalFit.C++(3,'${i}')' >log/PbPb3_pt${PTBIN[i]}to${PTBIN[i+1]}.log
# root -l -b -q 'SignalFit_wideRange.C++(3,'${i}')' 
root -l -b -q 'SignalFit_fixShapeCutScanTest.C++(3,'${i}')' >>temp.log 
i=$(($i+1))
done
