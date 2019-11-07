#!/bin/sh

#for i in {1..20}
#	do
#	ls /scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_Data/MB${i}/*.root > Ds_pp_Data_MB${i}.lis
#done

for j in {1..11}
	do
	ls /scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_Data/TrkOnly/HIMB${j}/*.root > Ds_PbPb_Data_HIMB${j}_TrkOnly.lis
done

ls /scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_Data/GJ/HIMB2_GJpart1/*.root  >Ds_PbPb_Data_HIMB2_GJ.lis
ls /scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_Data/GJ/HIMB2_GJpart2/*.root >>Ds_PbPb_Data_HIMB2_GJ.lis
ls /scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_Data/GJ/HIMB2_GJpart3/*.root >>Ds_PbPb_Data_HIMB2_GJ.lis
ls /scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_Data/GJ/HIMB2_GJpart4/*.root >>Ds_PbPb_Data_HIMB2_GJ.lis


ls /scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_Data/GJ/HIMB3_GJpart3/*.root  >Ds_PbPb_Data_HIMB3_GJ.lis
ls /scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_Data/GJ/HIMB3_GJpart4/*.root >>Ds_PbPb_Data_HIMB3_GJ.lis

ls /scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_Data/GJ/HIMB4_GJpart3/*.root  >Ds_PbPb_Data_HIMB4_GJ.lis
ls /scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_Data/GJ/HIMB4_GJpart4/*.root >>Ds_PbPb_Data_HIMB4_GJ.lis
