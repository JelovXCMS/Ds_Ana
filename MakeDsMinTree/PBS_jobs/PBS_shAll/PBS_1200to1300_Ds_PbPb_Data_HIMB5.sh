#!/bin/sh

#PBS -l nodes=1,walltime=04:00:00,mem=2400MB
#PBS -N Dnt1200to1300_Ds_PbPb_Data_HIMB5
#PBS -o /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/PBS_log/PBS_1200to1300_Ds_PbPb_Data_HIMB5.log
#PBS -e /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/PBS_log/PBS_1200to1300_Ds_PbPb_Data_HIMB5.err
#PBS -r n
#PBS -V
#PBS -q standby 

cd /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree

./exec_condorfilelist.sh input_list/Ds_PbPb_Data_HIMB5.lis 1 1200 1300 /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/cent30100/Ds_PbPb_Data_HIMB5 1

