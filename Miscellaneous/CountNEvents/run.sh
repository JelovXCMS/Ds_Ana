#!/bin/sh


for mbi in {1..11}
do
	root -b -q -l CountNEvents.C++\(\"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/input_list/Ds_PbPb_Data_HIMB${mbi}_TrkOnly.lis\"\) >log/PbPb_Data_HIMB${mbi}_TrkOnly.txt

done


for mbi2 in {2..4}
do
	root -b -q -l CountNEvents.C++\(\"/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/input_list/Ds_PbPb_Data_HIMB${mbi2}_GJ.lis\"\) >log/PbPb_Data_HIMB${mbi2}_GJ.txt

done
