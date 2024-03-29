#!/bin/sh

#PBS -l nodes=1,walltime=04:00:00,mem=3600MB
#PBS -N Dnt540to570_Ds_pp_Data_MB5
#PBS -o /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/PBS_log/PBS_540to570_Ds_pp_Data_MB5.log
#PBS -e /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/PBS_log/PBS_540to570_Ds_pp_Data_MB5.err
#PBS -r n
#PBS -V
#PBS -q standby 

cd /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree

./exec_condorfilelist.sh input_list/Ds_pp_Data_MB5_kpkptest.lis 0 540 570 /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data_samesign/Ds_pp_Data_MB5_kpkptest 1

