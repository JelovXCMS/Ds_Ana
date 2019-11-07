#!/bin/sh

#PBS -l nodes=1,walltime=04:00:00,mem=3600MB
#PBS -N Dnt40to80_Ds_pp_Data_MB14
#PBS -o /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/PBS_log/PBS_40to80_Ds_pp_Data_MB14.log
#PBS -e /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/PBS_log/PBS_40to80_Ds_pp_Data_MB14.err
#PBS -r n
#PBS -V
#PBS -q standby 

cd /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree

./exec_condorfilelist.sh input_list/Ds_pp_Data_MB14_kpkptest.lis 0 40 80 /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/Ds_pp_Data_MB14_kpkptest 1

