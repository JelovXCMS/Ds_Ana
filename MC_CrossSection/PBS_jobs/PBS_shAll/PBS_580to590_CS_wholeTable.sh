#!/bin/sh

#PBS -l nodes=1,walltime=04:00:00,mem=3600MB
#PBS -N Dnt580to590_CS_wholeTable
#PBS -o /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_CrossSection/PBS_jobs/PBS_log/PBS_580to590_CS_wholeTable.log
#PBS -e /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_CrossSection/PBS_jobs/PBS_log/PBS_580to590_CS_wholeTable.err
#PBS -r n
#PBS -V
#PBS -q standby 

cd /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_CrossSection

./exec_condorfilelist.sh input_list/CMStune_CharmAll.lis 0 580 590 /scratch/halstead/p/peng43/Ds_phikkpi/CS_CharmAll 1

