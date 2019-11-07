#!/bin/sh


#root -b -q BuildFitfile.C++
#root -b -q 'BuildMCfile.C++(0,0,0)'
#root -b -q 'BuildMCfile.C++(0,0,1)'
#root -b -q 'BuildMCfile.C++(0,1,0)'
#root -b -q 'BuildMCfile.C++(0,1,1)'

root -b -q 'Fit_sideband.C++ (0,"Ddls")'
root -b -q 'Fit_sideband.C++(0,"Dalpha")'
root -b -q 'Fit_sideband.C++(0,"Dchi2cl")'
