#!/bin/sh

#PBS -l nodes=1,walltime=04:00:00,mem=3600MB
#PBS -N Dnt3750to3900_Ds_PbPb_Data_HIMB3_GJ
#PBS -o /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/PBS_log/PBS_3750to3900_Ds_PbPb_Data_HIMB3_GJ.log
#PBS -e /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/PBS_log/PBS_3750to3900_Ds_PbPb_Data_HIMB3_GJ.err
#PBS -r n
#PBS -V
#PBS -q standby 

cd /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree

./exec_condorfilelist.sh input_list/Ds_PbPb_Data_HIMB3_GJ.lis 3 3750 3900 /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/GJ/Ds_PbPb_Data_HIMB3_GJ 1

