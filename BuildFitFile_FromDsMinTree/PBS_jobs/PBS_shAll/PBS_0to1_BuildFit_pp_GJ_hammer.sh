#!/bin/sh

#PBS -l nodes=1,walltime=04:00:00,mem=4800MB
#PBS -N Dnt0to1_BuildFit_pp_GJ_hammer
#PBS -o /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BuildFitFile_FromDsMinTree/PBS_jobs/PBS_log/PBS_0to1_BuildFit_pp_GJ_hammer.log
#PBS -e /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BuildFitFile_FromDsMinTree/PBS_jobs/PBS_log/PBS_0to1_BuildFit_pp_GJ_hammer.err
#PBS -r n
#PBS -V
#PBS -q standby 

cd /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BuildFitFile_FromDsMinTree

./exec_condorfilelist.sh input_list/pp_Data_DsMinTree.lis 0 0 1 /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BuildFitFile_FromDsMinTree/output_temp 1

