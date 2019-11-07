#!/bin/sh

#PBS -l nodes=1,walltime=04:00:00,mem=3600MB
#PBS -N Dnt470to480_CS_wholeTable_CharmAllph20
#PBS -o /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_CrossSection/PBS_jobs/PBS_log/PBS_470to480_CS_wholeTable_CharmAllph20.log
#PBS -e /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_CrossSection/PBS_jobs/PBS_log/PBS_470to480_CS_wholeTable_CharmAllph20.err
#PBS -r n
#PBS -V
#PBS -q standby 

cd /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_CrossSection

./exec_condorfilelist.sh input_list/CMStune_CharmAll_pthat20.lis 0 470 480 /scratch/halstead/p/peng43/Ds_phikkpi/CS_CharmAll_ph20 1

