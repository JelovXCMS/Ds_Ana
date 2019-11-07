#!/bin/sh

#PBS -l nodes=1,walltime=04:00:00,mem=3600MB
#PBS -N Dnt3600to3715_Ds_PbPb_Data_HIMB3_TrkOnly
#PBS -o /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/PBS_log/PBS_3600to3715_Ds_PbPb_Data_HIMB3_TrkOnly.log
#PBS -e /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/PBS_log/PBS_3600to3715_Ds_PbPb_Data_HIMB3_TrkOnly.err
#PBS -r n
#PBS -V
#PBS -q standby 

cd /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree

./exec_condorfilelist.sh input_list/Ds_PbPb_Data_HIMB3_TrkOnly.lis 3 3600 3715 /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/TrkOnly/Ds_PbPb_Data_HIMB3_TrkOnly 1

