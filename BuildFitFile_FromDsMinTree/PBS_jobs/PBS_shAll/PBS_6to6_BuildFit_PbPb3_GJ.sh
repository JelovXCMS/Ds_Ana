#!/bin/sh

#PBS -l nodes=1,walltime=04:00:00,mem=2400MB
#PBS -N Dnt6to6_BuildFit_PbPb3_GJ
#PBS -o /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BuildFitFile_FromDsMinTree/PBS_jobs/PBS_log/PBS_6to6_BuildFit_PbPb3_GJ.log
#PBS -e /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BuildFitFile_FromDsMinTree/PBS_jobs/PBS_log/PBS_6to6_BuildFit_PbPb3_GJ.err
#PBS -r n
#PBS -V
#PBS -q standby 

cd /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BuildFitFile_FromDsMinTree

./exec_condorfilelist.sh input_list/PbPb_Data_DsMinTree_MB2_GJ.lis 3 6 6 /scratch/halstead/p/peng43/Ds_phikkpi/BuildFitFile_FromDsMinTree/PbPb3_Data/ 1

