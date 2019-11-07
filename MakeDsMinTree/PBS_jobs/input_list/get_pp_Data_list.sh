#!/bin/sh

for i in {1..20}
do
#  ls /scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_Data/MB${i}/*.root >Ds_pp_Data_MB${i}.lis
  ls /mnt/hadoop/store/group/hi/chengchi/Dntuple/Dsntuple/pp_Data/MB${i}/*.root >Ds_pp_Data_MB${i}.lis
done

sed -i 's/\/mnt\/hadoop/root:\/\/xrootd.rcac.purdue.edu\//g' Ds_pp_Data_MB*.lis
