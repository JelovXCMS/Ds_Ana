#!/bin/sh

#PBS -l nodes=1,walltime=04:00:00,mem=3600MB
#PBS -N Dnt100to110_Ds_PbPb_MC_Prompt_Phi_pt6
#PBS -o /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/PBS_log/PBS_100to110_Ds_PbPb_MC_Prompt_Phi_pt6.log
#PBS -e /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree/PBS_jobs/PBS_log/PBS_100to110_Ds_PbPb_MC_Prompt_Phi_pt6.err
#PBS -r n
#PBS -V
#PBS -q standby 

cd /home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MakeDsMinTree

./exec_condorfilelist.sh input_list/Ds_PbPb_MC_Prompt_Phi_pt6.lis 1 100 110 /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_offMC/Ds_PbPb_MC_Prompt_Phi_pt6 0

