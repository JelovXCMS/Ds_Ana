#!/bin/sh

#PBS -l nodes=1,walltime=04:00:00,mem=3600MB
#PBS -N Dnt600to629_Ds_pp_Data_MB17
#PBS -o /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/PBS_log/PBS_600to629_Ds_pp_Data_MB17.log
#PBS -e /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/PBS_log/PBS_600to629_Ds_pp_Data_MB17.err
#PBS -r n
#PBS -V
#PBS -q standby 

cd /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree

./exec_condorfilelist.sh input_list/Ds_pp_Data_MB17.lis 0 600 629 /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/Ds_pp_Data_MB17 1

