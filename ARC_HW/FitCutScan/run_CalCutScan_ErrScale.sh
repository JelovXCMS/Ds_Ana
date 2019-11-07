#!/bin/bash


Runpp=1
RunPbPb=0

RunDdls=0
RunDalpha=1
RunDchi2cl=1
RunDphimass=0

rm CalCutScan_AnaBin_Smear.exe
g++ CalCutScan_AnaBin_Smear.C $(root-config --cflags --libs) -L $ROOTSYS/lib -lRooFit -lHtml -lMinuit -lRooFitCore -lRooStats   -Wall -O2 -o CalCutScan_AnaBin_Smear.exe

#int CalCutScan_AnaBin_Smear(int isPbPb=0, int useAnaBin=1, int ibin_Dpt=7, double DptLow=10, double DptHigh=20, TString var_scan="Dalpha", int FixShape=1, int useReweight=1, TString ReweightTree="DdxyzErrWeight", int doPFrScan=1, int fixPFr=0, double PfrVal=0.7){
ppPbPb=0
useAnaBin=1
ibin_Dpt_Low=0
ibin_Dpt_High=7
DptLow=2
DptHigh=3
FixShape=1
useReweight=0
ReweightTree=DdxyzErrWeight
doPFrScan=0
fixPFr=0
PfrVal=0.85

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

    if [ ${RunDdls} == 1 ]
    then
      echo "do RunDdls"
      ./CalCutScan_AnaBin_Smear.exe ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} Ddls ${FixShape} 0  ${ReweightTree} ${doPFrScan} ${fixPFr} ${PfrVal}
      ./CalCutScan_AnaBin_Smear.exe ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} Ddls ${FixShape} ${useReweight} ${ReweightTree} ${doPFrScan} ${fixPFr} ${PfrVal}
    fi
    if [ ${RunDalpha} == 1 ]
    then
      echo "do RunDalpha"
      ./CalCutScan_AnaBin_Smear.exe ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} Dalpha ${FixShape} 0  ${ReweightTree} ${doPFrScan} ${fixPFr} ${PfrVal}
      ./CalCutScan_AnaBin_Smear.exe ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} Dalpha ${FixShape} ${useReweight} ${ReweightTree} ${doPFrScan} ${fixPFr} ${PfrVal}
    fi
    if [ ${RunDchi2cl} == 1 ]
    then
      echo "do RunDchi2cl"
      ./CalCutScan_AnaBin_Smear.exe ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} Dchi2cl ${FixShape} 0  ${ReweightTree} ${doPFrScan} ${fixPFr} ${PfrVal}
      ./CalCutScan_AnaBin_Smear.exe ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} Dchi2cl ${FixShape} ${useReweight} ${ReweightTree} ${doPFrScan} ${fixPFr} ${PfrVal}
    fi
    if [ ${RunDphimass} == 1 ]
    then
      echo "do RunDphimass"
      ./CalCutScan_AnaBin_Smear.exe ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} PhiMass ${FixShape} 0  ${ReweightTree} ${doPFrScan} ${fixPFr} ${PfrVal}
      ./CalCutScan_AnaBin_Smear.exe ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} PhiMass ${FixShape} ${useReweight} ${ReweightTree} ${doPFrScan} ${fixPFr} ${PfrVal}
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

  echo "Runpp"
  while (( ${ibin_Dpt} < $((ibin_Dpt_High+1)) ))
  do
    echo "ibin_Dpt = "${ibin_Dpt}

    if [ ${RunDdls} == 1 ]
    then
      echo "do RunDdls"
      ./CalCutScan_AnaBin_Smear.exe ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} Ddls ${FixShape} 0  ${ReweightTree} ${doPFrScan} ${fixPFr} ${PfrVal}
      ./CalCutScan_AnaBin_Smear.exe ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} Ddls ${FixShape} ${useReweight} ${ReweightTree} ${doPFrScan} ${fixPFr} ${PfrVal}
    fi
    if [ ${RunDalpha} == 1 ]
    then
      echo "do RunDalpha"
      ./CalCutScan_AnaBin_Smear.exe ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} Dalpha ${FixShape} 0  ${ReweightTree} ${doPFrScan} ${fixPFr} ${PfrVal}
      ./CalCutScan_AnaBin_Smear.exe ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} Dalpha ${FixShape} ${useReweight} ${ReweightTree} ${doPFrScan} ${fixPFr} ${PfrVal}
    fi
    if [ ${RunDchi2cl} == 1 ]
    then
      echo "do RunDchi2cl"
      ./CalCutScan_AnaBin_Smear.exe ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} Dchi2cl ${FixShape} 0  ${ReweightTree} ${doPFrScan} ${fixPFr} ${PfrVal}
      ./CalCutScan_AnaBin_Smear.exe ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} Dchi2cl ${FixShape} ${useReweight} ${ReweightTree} ${doPFrScan} ${fixPFr} ${PfrVal}
    fi
    if [ ${RunDphimass} == 1 ]
    then
      echo "do RunDphimass"
      ./CalCutScan_AnaBin_Smear.exe ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} PhiMass ${FixShape} 0  ${ReweightTree} ${doPFrScan} ${fixPFr} ${PfrVal}
      ./CalCutScan_AnaBin_Smear.exe ${ppPbPb} ${useAnaBin} ${ibin_Dpt} ${DptLow} ${DptHigh} PhiMass ${FixShape} ${useReweight} ${ReweightTree} ${doPFrScan} ${fixPFr} ${PfrVal}
    fi

    ibin_Dpt=$((ibin_Dpt+1))
  done
  # ./run_Fit_AnaBin.sh 0 1 7 3 4 Ddls 1 1

fi


