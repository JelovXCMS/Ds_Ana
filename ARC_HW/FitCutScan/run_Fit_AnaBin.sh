#!/bin/bash

doDefaultFit=$9

isPbPb=$1
useAnaBin=$2
ibin_Dpt=$3
DptLow=$4
DptHigh=$5
# var_scan="Ddls"
var_scan=$6
# scan_val=$7
FixShape=$7
FixShapeVal=$8

scan_val=(1.5 1.75 2.0 2.25 2.5 2.75 3.0 3.25 3.5 3.75 4.0)
# scan_val=(3.5 4.0 4.5 5.0)
# scan_val=(3.25 3.5 3.75 4.0 4.25 4.5 4.75 5.0 5.25 5.5)
# scan_val=(2.75 3.0 3.25 3.5 3.75 4.0 4.25 4.5 4.75 5.0 5.25 5.5)
# scan_val=(2.5 2.75 3.0)
# scan_val=(2.0 3.5)

# DdlsScan_val=(1.5 1.75 2.0 2.25 2.5 2.75 3.0 3.25 3.5 3.75 4.0)
# DdlsScan_val=(1.5 2.5 3.5 4.5 5.5)

# use default scan array in code if useAnaBin, ignore the following

DdlsScan_pp_val=(1.5 2.5 3.5 4.5 5.5) 
DdlsScan_PbPb_val=(2.5 3.5 4.5 5.5 6.5) 
DdlsScan_PbPb_Low_val=(3.5 4.5 5.5 6.5) 

Dchi2clScan_pp_val=(0.02 0.12 0.22 0.32 0.42)
Dchi2clScan_PbPb_val=(0.05 0.15 0.25 0.35 0.45)

DalphaScan_val=(0.08 0.11 0.14 0.17 0.2)

PhiMassScan_val=(0.008 0.0085 0.009 0.0095 0.01 0.0105 0.0110)

jobscount=0

######
#if run default fit
#then
#run default fit
#
######

#int FitCutScan_AnaBin(int isPbPb=0,int useAnaBin=1, int ibin_Dpt=7, double DptLow=6, double DptHigh=40, TString var_scan="Dchi2cl", double scan_val=0.05, int FixShape=0, double shapeVal=3.0 , int useDefaultCut=0){


#rm FitCutScan_AnaBin.exe
#g++ FitCutScan_AnaBin.C $(root-config --cflags --libs) -L $ROOTSYS/lib -lRooFit -lHtml -lMinuit -lRooFitCore -lRooStats   -Wall -O2 -o FitCutScan_AnaBin.exe

if [ ${doDefaultFit} -eq 1 ]
then
	echo "do default Fit"
	./FitCutScan_AnaBin.exe  ${isPbPb} ${useAnaBin} ${ibin_Dpt}  ${DptLow}  ${DptHigh}  ${var_scan} 0 0 ${FixShapeVal} 1  >log/Fit_isPbPb${isPbPb}_ibinDpt${ibin_Dpt}_DefaultCutFit.log &
	wait
	./FitCutScan_AnaBin.exe  ${isPbPb} ${useAnaBin} ${ibin_Dpt}  ${DptLow}  ${DptHigh}  ${var_scan} 0 1 ${FixShapeVal} 1  >log/Fit_isPbPb${isPbPb}_ibinDpt${ibin_Dpt}_DefaultCutFit.log &
	wait 
fi



if [ ${var_scan} == 'Ddls' ] && [ ${isPbPb} -eq 0 ] 
then 

	echo "Ddls pp "

for i in ${DdlsScan_pp_val[@]}
do
	echo ${var_scan} " = " ${i}

	./FitCutScan_AnaBin.exe  ${isPbPb} ${useAnaBin} ${ibin_Dpt}  ${DptLow}  ${DptHigh}  ${var_scan}  ${i}  ${FixShape} ${FixShapeVal} 0 >log/Fit_isPbPb${isPbPb}_ibinDpt${ibin_Dpt}_${var_scan}${i}_FixShape${FixShape}.log &

	 jobscount=$((${jobscount}+1))
	echo "jobscount = "${jobscount}
	if [ $((${jobscount} % 2)) -eq 0 ]
	then
		echo "wait"
		wait
	fi


done
fi # end if var=Ddls

wait



if [ ${var_scan} == 'Ddls' ] && [ ${isPbPb} -gt 0 ] && [ ${ibin_Dpt} -lt 3 ]
then 
		echo "PbPb low pt Ddls , ibin_Dpt = " ${ibin_Dpt}


for i in ${DdlsScan_PbPb_Low_val[@]}
do
	echo ${var_scan} " = " ${i}

	./FitCutScan_AnaBin.exe  ${isPbPb} ${useAnaBin} ${ibin_Dpt}  ${DptLow}  ${DptHigh}  ${var_scan}  ${i}  ${FixShape} ${FixShapeVal} 0 >log/Fit_isPbPb${isPbPb}_ibinDpt${ibin_Dpt}_${var_scan}${i}_FixShape${FixShape}.log &

	 jobscount=$((${jobscount}+1))
	echo "jobscount = "${jobscount}
	if [ $((${jobscount} % 2)) -eq 0 ]
	then
		echo "wait"
		wait
	fi


done
fi # end if var=Ddls

wait


if [ ${var_scan} == 'Ddls' ] && [ ${isPbPb} -gt 0 ] && [ ${ibin_Dpt} -gt 2 ]
then 
		echo "PbPb high pt Ddls , ibin_Dpt = " ${ibin_Dpt}


for i in ${DdlsScan_PbPb_val[@]}
do
	echo ${var_scan} " = " ${i}

	./FitCutScan_AnaBin.exe  ${isPbPb} ${useAnaBin} ${ibin_Dpt}  ${DptLow}  ${DptHigh}  ${var_scan}  ${i}  ${FixShape} ${FixShapeVal} 0 >log/Fit_isPbPb${isPbPb}_ibinDpt${ibin_Dpt}_${var_scan}${i}_FixShape${FixShape}.log &

	 jobscount=$((${jobscount}+1))
	echo "jobscount = "${jobscount}
	if [ $((${jobscount} % 2)) -eq 0 ]
	then
		echo "wait"
		wait
	fi


done
fi # end if var=Ddls

wait





if [ ${var_scan} == 'Dalpha' ]
then 

for i in ${DalphaScan_val[@]}
do
	echo ${var_scan} " = " ${i}

	./FitCutScan_AnaBin.exe  ${isPbPb} ${useAnaBin} ${ibin_Dpt}  ${DptLow}  ${DptHigh}  ${var_scan}  ${i}  ${FixShape} ${FixShapeVal} 0 >log/Fit_isPbPb${isPbPb}_ibinDpt${ibin_Dpt}_${var_scan}${i}_FixShape${FixShape}.log &

	 jobscount=$((${jobscount}+1))
	echo "jobscount = "${jobscount}
	if [ $((${jobscount} % 2)) -eq 0 ]
	then
		echo "wait"
		wait
	fi


done
fi # end if var=Ddls

wait


if [ ${var_scan} == 'Dchi2cl' ] && [ ${isPbPb} -eq 0 ] 
then 

for i in ${Dchi2clScan_pp_val[@]}
do
	echo ${var_scan} " = " ${i}

	./FitCutScan_AnaBin.exe  ${isPbPb} ${useAnaBin} ${ibin_Dpt}  ${DptLow}  ${DptHigh}  ${var_scan}  ${i}  ${FixShape} ${FixShapeVal} 0 >log/Fit_isPbPb${isPbPb}_ibinDpt${ibin_Dpt}_${var_scan}${i}_FixShape${FixShape}.log &

	 jobscount=$((${jobscount}+1))
	echo "jobscount = "${jobscount}
	if [ $((${jobscount} % 2)) -eq 0 ]
	then
		echo "wait"
		wait
	fi


done
fi # end if var=Ddls

wait



if [ ${var_scan} == 'Dchi2cl' ] && [ ${isPbPb} -gt 0 ] 
then 

for i in ${Dchi2clScan_PbPb_val[@]}
do
	echo ${var_scan} " = " ${i}

	./FitCutScan_AnaBin.exe  ${isPbPb} ${useAnaBin} ${ibin_Dpt}  ${DptLow}  ${DptHigh}  ${var_scan}  ${i}  ${FixShape} ${FixShapeVal} 0 >log/Fit_isPbPb${isPbPb}_ibinDpt${ibin_Dpt}_${var_scan}${i}_FixShape${FixShape}.log &

	 jobscount=$((${jobscount}+1))
	echo "jobscount = "${jobscount}
	if [ $((${jobscount} % 2)) -eq 0 ]
	then
		echo "wait"
		wait
	fi


done
fi # end if var=Ddls

wait




if [ ${var_scan} == 'PhiMass' ]
then 

for i in ${PhiMassScan_val[@]}
do
	echo ${var_scan} " = " ${i}

	./FitCutScan_AnaBin.exe  ${isPbPb} ${useAnaBin} ${ibin_Dpt}  ${DptLow}  ${DptHigh}  ${var_scan}  ${i}  ${FixShape} ${FixShapeVal} 0 >log/Fit_isPbPb${isPbPb}_ibinDpt${ibin_Dpt}_${var_scan}${i}_FixShape${FixShape}.log &

	 jobscount=$((${jobscount}+1))
	echo "jobscount = "${jobscount}
	if [ $((${jobscount} % 2)) -eq 0 ]
	then
		echo "wait"
		wait
	fi


done
fi # end if var=Ddls

wait




