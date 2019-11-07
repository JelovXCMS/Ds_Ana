#!/bin/sh

#PBS -l nodes=1,walltime=04:00:00,mem=3600MB
#PBS -N Dnt0to100_Ds_PbPb_Data_HIMB2_TrkOnly
#PBS -o /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/PBS_log/PBS_0to100_Ds_PbPb_Data_HIMB2_TrkOnly.log
#PBS -e /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/PBS_log/PBS_0to100_Ds_PbPb_Data_HIMB2_TrkOnly.err
#PBS -r n
#PBS -V
#PBS -q standby 

cd /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree

./exec_condorfilelist.sh input_list/Ds_PbPb_Data_HIMB2_TrkOnly.lis 3 0 100 /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/TrkOnly/Ds_PbPb_Data_HIMB2_TrkOnly 1

