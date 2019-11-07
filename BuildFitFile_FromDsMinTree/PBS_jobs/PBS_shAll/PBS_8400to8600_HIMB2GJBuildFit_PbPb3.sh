#!/bin/sh

#PBS -l nodes=1,walltime=04:00:00,mem=4800MB,naccesspolicy=singleuser
#PBS -N Dnt8400to8600_HIMB2GJBuildFit_PbPb3
#PBS -o /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BuildFitFile_FromDsMinTree/PBS_jobs/PBS_log/PBS_8400to8600_HIMB2GJBuildFit_PbPb3.log
#PBS -e /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BuildFitFile_FromDsMinTree/PBS_jobs/PBS_log/PBS_8400to8600_HIMB2GJBuildFit_PbPb3.err
#PBS -r n
#PBS -V
#PBS -q standby 

cd /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BuildFitFile_FromDsMinTree

./exec_condorfilelist.sh input_list/Ds_PbPb_Data_HIMB2_GJ.lis 3 8400 8600 /scratch/halstead/p/peng43/Ds_phikkpi/BuildFitFile_FromDsMinTree/PbPb3_Data_moreScan//HIMB2_GJ/ 1
