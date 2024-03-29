#!/bin/sh

#PBS -l nodes=1,walltime=04:00:00,mem=4800MB,naccesspolicy=singleuser
#PBS -N Dnt6000to6200_HIMB4GJBuildFit_PbPb3
#PBS -o /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BuildFitFile_FromDsMinTree/PBS_jobs/PBS_log/PBS_6000to6200_HIMB4GJBuildFit_PbPb3.log
#PBS -e /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BuildFitFile_FromDsMinTree/PBS_jobs/PBS_log/PBS_6000to6200_HIMB4GJBuildFit_PbPb3.err
#PBS -r n
#PBS -V
#PBS -q standby 

cd /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BuildFitFile_FromDsMinTree

./exec_condorfilelist.sh input_list/Ds_PbPb_Data_HIMB4_GJ.lis 3 6000 6200 /scratch/halstead/p/peng43/Ds_phikkpi/BuildFitFile_FromDsMinTree/PbPb3_Data_moreScan//HIMB4_GJ/ 1

