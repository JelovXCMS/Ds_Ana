#!/bin/bash


cd /home/peng43/root_from_CMSSW/root610_09/CMSSW_10_1_0_pre1/src
eval `scramv1 runtime -sh`
cd -

rm Fit_sideband.exe
g++ Fit_sideband.C $(root-config --cflags --libs) -L $ROOTSYS/lib -lRooFit -lHtml -lMinuit -lRooFitCore -Wall -O2 -o Fit_sideband.exe

rm plot_compare.exe
g++ plot_compare.C $(root-config --cflags --libs)  -Wall -O2 -o plot_compare.exe

isPbPb=3
var_compare="DdxyzErr"
var_cut=("Dtrk1Pt" "Dtrk2Pt" "Dtrk3Pt")
var_cutLow=0.75
var_cutHigh=1.0


nvar_cut=${#var_cut[@]}

echo "var_cut 1 = " ${var_cut[1]}

# var_cutbins=(0.75 1.0 1.25 1.5 1.75 2 2.25 2.5 2.75)
# var_cutbins=(0.75 2 3 4 5 6 7 8 10 12 15)
var_cutbins=(1 2 3 4 5 6 8 12)
# var_cutbins=(0.75 2 3 4 5 7 10 15)
# var_cutbins=(0.75 1.0 )


var_bins=${#var_cutbins[@]}

echo ${var_bins}

mkdir -p log

# source compile.sh

# for i in {0..${var_bins} }
for (( j=0; j<nvar_cut; j++  ))
do
	echo "do varcut " ${var_cut[${j}]} 
	for (( i=0; i<${var_bins}-1; i++   ))
# for (( i=0; i<5 ; i++   ))
	do
		echo ${var_cutbins[($i)]} "to" ${var_cutbins[($i+1)]}
 	./Fit_sideband.exe ${isPbPb} ${var_compare} ${var_cut[${j}]} ${var_cutbins[($i)]} ${var_cutbins[($i+1)]} > log/Fit_sideband_pPb${isPbPb}_${var_compare}_${var_cut[${j}]}_${var_cutbins[($i)]}_${var_cutbins[($i+1)]}.log
	done
		./plot_compare.exe ${isPbPb} ${var_compare} ${var_cut[${j}]}
done

# ./plot_compare.exe ${isPbPb} ${var_compare} ${var_cut}

# echo ${#var_cutbins[@]}
# for i in 

# ./Fit_sideband.exe ${isPbPb} ${var_compare} ${var_cut} ${var_cutLow} ${var_cutHigh}

