#!/bin/sh

#PBS -l nodes=1,walltime=04:00:00,mem=2400MB
#PBS -N Dnt400to423_Ds_pp_Data_MB2
#PBS -o /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/PBS_log/PBS_400to423_Ds_pp_Data_MB2.log
#PBS -e /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/PBS_log/PBS_400to423_Ds_pp_Data_MB2.err
#PBS -r n
#PBS -V
#PBS -q standby 

cd /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree

./exec_condorfilelist.sh input_list/Ds_pp_Data_MB2.lis 0 400 423 /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB2 1

