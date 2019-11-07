#!/bin/sh

#PBS -l nodes=1,walltime=04:00:00,mem=2400MB
#PBS -N Dnt690to700_Ds_PbPb_Data_HIMB11_TrkOnly
#PBS -o /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/PBS_log/PBS_690to700_Ds_PbPb_Data_HIMB11_TrkOnly_V2.log
#PBS -e /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/PBS_log/PBS_690to700_Ds_PbPb_Data_HIMB11_TrkOnly_V2.err
#PBS -r n
#PBS -V
#PBS -q standby 

cd /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree

./exec_condorfilelist_V2.sh input_list/Ds_PbPb_Data_HIMB11_TrkOnly.lis 3 690 700 /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTreeV2/PbPb_Data/TrkOnly/Ds_PbPb_Data_HIMB11_TrkOnly 1

