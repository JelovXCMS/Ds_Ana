#!/bin/sh

#PBS -l nodes=1,walltime=04:00:00,mem=3600MB
#PBS -N Dnt3280to3290_CS_wholeTable_softQCD
#PBS -o /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_CrossSection/PBS_jobs/PBS_log/PBS_3280to3290_CS_wholeTable_softQCD.log
#PBS -e /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_CrossSection/PBS_jobs/PBS_log/PBS_3280to3290_CS_wholeTable_softQCD.err
#PBS -r n
#PBS -V
#PBS -q standby 

cd /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_CrossSection

./exec_condorfilelist.sh input_list/CMStune_CharmAll_softQCD.lis 0 3280 3290 /scratch/halstead/p/peng43/Ds_phikkpi/CS_CharmAll_softQCD 1

