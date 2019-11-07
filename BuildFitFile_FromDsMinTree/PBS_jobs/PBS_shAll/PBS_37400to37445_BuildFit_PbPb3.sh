#!/bin/sh

#PBS -l nodes=1,walltime=04:00:00,mem=4800MB,naccesspolicy=singleuser
#PBS -N Dnt37400to37445_BuildFit_PbPb3
#PBS -o /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BuildFitFile_FromDsMinTree/PBS_jobs/PBS_log/PBS_37400to37445_BuildFit_PbPb3.log
#PBS -e /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BuildFitFile_FromDsMinTree/PBS_jobs/PBS_log/PBS_37400to37445_BuildFit_PbPb3.err
#PBS -r n
#PBS -V
#PBS -q standby 

cd /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BuildFitFile_FromDsMinTree

./exec_condorfilelist.sh input_list/Ds_PbPb_Data_HIMB_GJ.lis 3 37400 37445 /scratch/halstead/p/peng43/Ds_phikkpi/BuildFitFile_FromDsMinTree/PbPb3_Data/ 1
