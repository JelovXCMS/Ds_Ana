#!/bin/bash


#isPbPb=$1
#useAnaBin=$2
#ibin_Dpt=$3
#DptLow=$4
#DptHigh=$5
## var_scan="Ddls"
#var_scan=$6
#FixShape=$7
#FixShapeVal=$8

Runpp=1
RunPbPb=1

RunDdls=0
RunDalpha=0
RunDchi2cl=0
RunDphimass=1

#  ./FitCutScan_AnaBin.exe  ${isPbPb} ${useAnaBin} ${ibin_Dpt}  ${DptLow}  ${DptHigh}  ${var_scan} ${i} ${FixShape} ${FixShapeVal} 1  >log/Fit_isPbPb${isPbPb}_ibinDpt${ibin_Dpt}_DefaultCutFit.log &
#isPbPb=$1
#useAnaBin=$2
#ibin_Dpt=$3
#DptLow=$4
#DptHigh=$5
#var_scan=$6
#FixShape=$7
#FixShapeVal=$8
#doDefaultFit=$9

useAnaBin=1
DptLow=2
DptHigh=4
ibin_Dpt_Low=0  #not used
ibin_Dpt_High=7 #not used
var_scan=Ddls
FixShape=1
FixShapeVal=1 # not used

doDefaultFit=1

ibin_Dpt=${ibin_Dpt_Low}

ppPbPb=0
# ./run_Fit_AnaBin.sh 0 1 6 3 4 Ddls 1 1 
# ./run_Fit_AnaBin.sh 0 1 7 3 4 Ddls 1 1 1 


rm FitCutScan_AnaBin.exe
g++ FitCutScan_AnaBin.C $(root-config --cflags --libs) -L $ROOTSYS/lib -lRooFit -lHtml -lMinuit -lRooFitCore -lRooStats   -Wall -O2 -o FitCutScan_AnaBin.exe


# ./run_Fit_AnaBin.sh 0 1 7 3 4 Dalpha 1 1 1 

if [ ${Runpp} == 1 ]
then
	ppPbPb=0
	if ((ibin_Dpt_Low < 0 )) 
	then 
		ibin_Dpt_Low=0
	fi
	if ((ibin_Dpt_High > 7 ))
	then 
		ibin_Dpt_High=7
	fi

	ibin_Dpt=${ibin_Dpt_Low}
	echo "ibin_Dpt_Low = " ${ibin_Dpt_Low} " , ibin_Dpt_High = " ${ibin_Dpt_High}


	echo "Runpp"
	while (( ${ibin_Dpt} < $((ibin_Dpt_High+1)) ))
	do 
		echo "ibin_Dpt = "${ibin_Dpt}
		doDefaultFit=1

		if [ ${RunDdls} == 1 ]
		then
			echo "do RunDdls"
			./run_Fit_AnaBin.sh ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} Ddls ${FixShape} ${FixShapeVal} ${doDefaultFit}
			doDefaultFit=0 
		fi
		if [ ${RunDalpha} == 1 ]
		then
			echo "do RunDalpha"
			./run_Fit_AnaBin.sh ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} Dalpha ${FixShape} ${FixShapeVal} ${doDefaultFit}
			doDefaultFit=0 
		fi
		if [ ${RunDchi2cl} == 1 ]
		then
			echo "do RunDchi2cl"
			./run_Fit_AnaBin.sh ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} Dchi2cl ${FixShape} ${FixShapeVal} ${doDefaultFit}
			doDefaultFit=0 
		fi
		if [ ${RunDphimass} == 1 ]
		then
			echo "do RunDphimass"
			./run_Fit_AnaBin.sh ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} PhiMass ${FixShape} ${FixShapeVal} ${doDefaultFit}
			doDefaultFit=0 
		fi

		ibin_Dpt=$((ibin_Dpt+1))
	done
	# ./run_Fit_AnaBin.sh 0 1 7 3 4 Ddls 1 1

fi


if [ ${RunPbPb} == 1 ]
then
	ppPbPb=3
	if ((ibin_Dpt_Low < 2 ))
	then 
		ibin_Dpt_Low=2
	fi
	if ((ibin_Dpt_High > 5 ))
	then 
		ibin_Dpt_High=5
	fi

	ibin_Dpt=${ibin_Dpt_Low}
	echo "ibin_Dpt_Low = " ${ibin_Dpt_Low} " , ibin_Dpt_High = " ${ibin_Dpt_High}


	echo "RunPbPb"
	while (( ${ibin_Dpt} < $((ibin_Dpt_High+1)) ))
	do 
		echo "ibin_Dpt = "${ibin_Dpt}
		doDefaultFit=1

		if [ ${RunDdls} == 1 ]
		then
			echo "do RunDdls"
			./run_Fit_AnaBin.sh ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} Ddls ${FixShape} ${FixShapeVal} ${doDefaultFit}
			doDefaultFit=0 
		fi
		if [ ${RunDalpha} == 1 ]
		then
			echo "do RunDalpha"
			./run_Fit_AnaBin.sh ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} Dalpha ${FixShape} ${FixShapeVal} ${doDefaultFit}
			doDefaultFit=0 
		fi
		if [ ${RunDchi2cl} == 1 ]
		then
			echo "do RunDchi2cl"
			./run_Fit_AnaBin.sh ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} Dchi2cl ${FixShape} ${FixShapeVal} ${doDefaultFit}
			doDefaultFit=0 
		fi
		if [ ${RunDphimass} == 1 ]
		then
			echo "do RunDphimass"
			./run_Fit_AnaBin.sh ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} PhiMass ${FixShape} ${FixShapeVal} ${doDefaultFit}
			doDefaultFit=0 
		fi

		ibin_Dpt=$((ibin_Dpt+1))
	done

fi











