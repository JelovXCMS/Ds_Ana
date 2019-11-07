#!/bin/bash

PbPb=3
PtLow=6
PtHigh=40
DVtxPCut=0.25
DdlsCut=4.5

rm ProjectToHis_PbPb.exe
g++ ProjectToHis_PbPb.C $(root-config --cflags --libs) -Wall -O2 -o ProjectToHis_PbPb.exe

# root -b -q ProjectToHis_PbPb.C++(${PbPb},${PtLow},${PtHigh},${DVtxPCut},${DdlsCut},"/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/GJ/Ds_PbPb_Data_HIMB2_GJ.root", "PbPb_MB2_GJ")
for i in {2..4}
do
	echo  "${PbPb} ${PtLow} ${PtHigh} ${DVtxPCut} ${DdlsCut}  /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/GJ/Ds_PbPb_Data_HIMB${i}_GJ.root"
	./ProjectToHis_PbPb.exe ${PbPb} ${PtLow} ${PtHigh} ${DVtxPCut} ${DdlsCut}  /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/GJ/Ds_PbPb_Data_HIMB${i}_GJ.root  PbPb_MB${i}_GJ  &
done

for i in {1..11}
do
	echo  "${PbPb} ${PtLow} ${PtHigh} ${DVtxPCut} ${DdlsCut}  /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/TrkOnly/Ds_PbPb_Data_HIMB${i}_TrkOnly.root"
	./ProjectToHis_PbPb.exe ${PbPb} ${PtLow} ${PtHigh} ${DVtxPCut} ${DdlsCut}  /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/TrkOnly/Ds_PbPb_Data_HIMB${i}_TrkOnly.root  PbPb_MB${i}_TrkOnly  &
done

wait
