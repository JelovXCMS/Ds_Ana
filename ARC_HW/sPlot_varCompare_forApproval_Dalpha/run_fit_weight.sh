#!/bin/bash


cd /home/peng43/root_from_CMSSW/root610_09/CMSSW_10_1_0_pre1/src
eval `scramv1 runtime -sh`
cd -

rm Fit_splot_onlySave.exe
# g++ Fit_splot_onlySave.C $(root-config --cflags --libs) -L $ROOTSYS/lib -lRooFit -lHtml -lMinuit -lRooFitCore -Wall -O2 -o Fit_splot_onlySave.exe
g++ Fit_splot_onlySave.C $(root-config --cflags --libs) -L $ROOTSYS/lib -lRooStats -lRooFit -lHtml -lMinuit -lRooFitCore -Wall -O2 -o Fit_splot_onlySave.exe
rm CompAndTune_DataMC.exe
g++ CompAndTune_DataMC.C $(root-config --cflags --libs)  -Wall -O2 -o CompAndTune_DataMC.exe

doFit=0
doCompWt=1

isPbPb=0

var_cutbins=(2 4 6 8 10 40)
#var_cutbins=(6 8 10 40)



var_compare="Dalpha"
var_cp_Low=0.0
var_cp_High=0.2
var_cp_DrawLow=0.0
var_cp_DrawHigh=0.2

#var_compare="Dchi2cl"
#var_cp_Low=0
#var_cp_High=1
#var_cp_DrawLow=0
#var_cp_DrawHigh=1

#var_compare="DdlErr"
#var_cp_Low=0
#var_cp_High=1
#var_cp_DrawLow=0
#var_cp_DrawHigh=0.1

#var_compare="Ddl"
#var_cp_Low=0
#var_cp_High=10
#var_cp_DrawLow=0
#var_cp_DrawHigh=1

#var_compare="Ddls"
#var_cp_Low=0
#var_cp_High=200
#var_cp_DrawLow=0
#var_cp_DrawHigh=30

# var_compare="DdxyzErr"
# var_cp_Low=0
# var_cp_High=1
# var_cp_DrawLow=0
# var_cp_DrawHigh=0.1



# var_cut=("Dtrk1Pt" "Dtrk2Pt" "Dtrk3Pt")
var_cut=("Dpt")

# var_cutLow=0.75
# var_cutHigh=1.0

nvar_cut=${#var_cut[@]}

echo "var_cut 1 = " ${var_cut[1]}

# var_cutbins=(0.75 1.0 1.25 1.5 1.75 2 2.25 2.5 2.75)
# var_cutbins=(0.75 2 3 4 5 6 7 8 10 12 15)
# var_cutbins=(6 10 40)
# var_cutbins=(6 8 10 20 40)
# var_cutbins=(6 8 10)
# var_cutbins=(6 8)
#var_cutbins=(2 4 6 8 10 40)
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
		if [ $doFit -eq 1 ]
			then 	
  			./Fit_splot_onlySave.exe ${isPbPb} ${var_compare} ${var_cut[${j}]} ${var_cutbins[($i)]} ${var_cutbins[($i+1)]} ${var_cp_Low} ${var_cp_High} ${var_cp_DrawLow} ${var_cp_DrawHigh} > log/Fit_splot_pPb${isPbPb}_${var_compare}_${var_cut[${j}]}_${var_cutbins[($i)]}_${var_cutbins[($i+1)]}_saveonly.log &
		fi
		# if [ $doCompWt -eq 1 ]
			# then
  		# ./CompAndTune_DataMC.exe ${isPbPb} ${var_compare} ${var_cut[${j}]} ${var_cutbins[($i)]} ${var_cutbins[($i+1)]} > log/CompAndTune_DataMC_pPb${isPbPb}_${var_compare}_${var_cut[${j}]}_${var_cutbins[($i)]}_${var_cutbins[($i+1)]}.log
		# fi
	done
		# ./plot_compare.exe ${isPbPb} ${var_compare} ${var_cut[${j}]}
done

wait
echo "fit done"

for (( j=0; j<nvar_cut; j++  ))
do
	echo "do varcut " ${var_cut[${j}]} 
	for (( i=0; i<${var_bins}-1; i++   ))
# for (( i=0; i<5 ; i++   ))
	do
		echo ${var_cutbins[($i)]} "to" ${var_cutbins[($i+1)]}
		##if [ $doFit -eq 1 ]
		##	then 	
  	##		./Fit_splot_onlySave.exe ${isPbPb} ${var_compare} ${var_cut[${j}]} ${var_cutbins[($i)]} ${var_cutbins[($i+1)]} ${var_cp_Low} ${var_cp_High} ${var_cp_DrawLow} ${var_cp_DrawHigh} > log/Fit_splot_pPb${isPbPb}_${var_compare}_${var_cut[${j}]}_${var_cutbins[($i)]}_${var_cutbins[($i+1)]}_saveonly.log &
		##fi
		if [ $doCompWt -eq 1 ]
			then
  		./CompAndTune_DataMC.exe ${isPbPb} ${var_compare} ${var_cut[${j}]} ${var_cutbins[($i)]} ${var_cutbins[($i+1)]} > log/CompAndTune_DataMC_pPb${isPbPb}_${var_compare}_${var_cut[${j}]}_${var_cutbins[($i)]}_${var_cutbins[($i+1)]}.log
		fi
	done
		# ./plot_compare.exe ${isPbPb} ${var_compare} ${var_cut[${j}]}
done

wait

	echo "Comp done"

# ./plot_compare.exe ${isPbPb} ${var_compare} ${var_cut}

# echo ${#var_cutbins[@]}
# for i in 

# ./Fit_sideband.exe ${isPbPb} ${var_compare} ${var_cut} ${var_cutLow} ${var_cutHigh}

