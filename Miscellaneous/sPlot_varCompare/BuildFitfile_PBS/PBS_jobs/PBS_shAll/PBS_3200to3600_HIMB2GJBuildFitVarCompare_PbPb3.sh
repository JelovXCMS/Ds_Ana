#!/bin/sh

#PBS -l nodes=1,walltime=04:00:00,mem=4800MB,naccesspolicy=singleuser
#PBS -N Dnt3200to3600_HIMB2GJBuildFitVarCompare_PbPb3
#PBS -o /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/Miscellaneous/sPlot_varCompare/BuildFitfile_PBS/PBS_jobs/PBS_log/PBS_3200to3600_HIMB2GJBuildFitVarCompare_PbPb3.log
#PBS -e /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/Miscellaneous/sPlot_varCompare/BuildFitfile_PBS/PBS_jobs/PBS_log/PBS_3200to3600_HIMB2GJBuildFitVarCompare_PbPb3.err
#PBS -r n
#PBS -V
#PBS -q standby 

cd /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/Miscellaneous/sPlot_varCompare/BuildFitfile_PBS

./exec_condorfilelist.sh input_list/Ds_PbPb_Data_HIMB2_GJ.lis 3 3200 3600 /scratch/halstead/p/peng43/Ds_phikkpi/BuildFitFile_ForVarCompare/PbPb3_Data//HIMB2_GJ/ 1

