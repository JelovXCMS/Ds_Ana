#!/bin/sh


#root -b -q BuildFitfile.C++
# must use PBS to generate

#root -b -q 'BuildMCfile.C++(3,0,0)'
#root -b -q 'BuildMCfile.C++(3,0,1)'
#root -b -q 'BuildMCfile.C++(3,1,0)'
#root -b -q 'BuildMCfile.C++(3,1,1)'

root -b -q 'Fit_sideband.C++ (3,"Ddls")'
root -b -q 'Fit_sideband.C++(3,"Dalpha")'
root -b -q 'Fit_sideband.C++(3,"Dchi2cl")'
