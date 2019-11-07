#!/bin/sh

for mbi in {1..20}
do
	./mutiCopy.sh /scratch/halstead/p/peng43/Ds_phikkpi/FileList/pp_Data_MB${mbi}.lis  Dntuple/pp_Data/MB${mbi} &>log/pp_Data_MB${mbi}.log   &
done

for mbi2 in {1..11}
do
	./mutiCopy.sh /scratch/halstead/p/peng43/Ds_phikkpi/FileList/PbPb_Data_TrkOnly_HIMB${mbi2}.lis  Dntuple/PbPb_Data/TrkOnly/HIMB${mbi2} &>log/PbPb_Data_TrkOnly_MB${mbi2}.log   &
done

#	./mutiCopy.sh /scratch/halstead/p/peng43/Ds_phikkpi/FileList/PbPb_Data_GJ_HIMB2_part1.lis Dntuple/PbPb_Data/GJ/HIMB2_GJpart1 &>log/PbPb_Data_GJ_HIMB2_part1.log
	./mutiCopy.sh /scratch/halstead/p/peng43/Ds_phikkpi/FileList/PbPb_Data_GJ_HIMB2_part2.lis Dntuple/PbPb_Data/GJ/HIMB2_GJpart2 &>log/PbPb_Data_GJ_HIMB2_part2.log
	./mutiCopy.sh /scratch/halstead/p/peng43/Ds_phikkpi/FileList/PbPb_Data_GJ_HIMB2_part3.lis Dntuple/PbPb_Data/GJ/HIMB2_GJpart3 &>log/PbPb_Data_GJ_HIMB2_part3.log
	./mutiCopy.sh /scratch/halstead/p/peng43/Ds_phikkpi/FileList/PbPb_Data_GJ_HIMB2_part4.lis Dntuple/PbPb_Data/GJ/HIMB2_GJpart4 &>log/PbPb_Data_GJ_HIMB2_part4.log


	./mutiCopy.sh /scratch/halstead/p/peng43/Ds_phikkpi/FileList/PbPb_Data_GJ_HIMB3_part3.lis Dntuple/PbPb_Data/GJ/HIMB3_GJpart3 &>log/PbPb_Data_GJ_HIMB3_part3.log
	./mutiCopy.sh /scratch/halstead/p/peng43/Ds_phikkpi/FileList/PbPb_Data_GJ_HIMB3_part4.lis Dntuple/PbPb_Data/GJ/HIMB3_GJpart4 &>log/PbPb_Data_GJ_HIMB3_part4.log

	./mutiCopy.sh /scratch/halstead/p/peng43/Ds_phikkpi/FileList/PbPb_Data_GJ_HIMB4_part3.lis Dntuple/PbPb_Data/GJ/HIMB4_GJpart3 &>log/PbPb_Data_GJ_HIMB4_part3.log
	./mutiCopy.sh /scratch/halstead/p/peng43/Ds_phikkpi/FileList/PbPb_Data_GJ_HIMB4_part4.lis Dntuple/PbPb_Data/GJ/HIMB4_GJpart4 &>log/PbPb_Data_GJ_HIMB4_part4.log
