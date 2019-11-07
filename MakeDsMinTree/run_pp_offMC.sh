#!/bin/bash

#                                     $input                  $output  isREAL  isPbPb
#root -b -l -q 'MakeDsMinTree.C++ ("Ds_phikkpi_pt4.root" ,"test.root", 0 ,0)'

mkdir -p /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_offMC/

ptbins=(0 2 4 6 10 19)
echo here

for i in {0..5}
do
#	echo $i
#	root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/offpp_MC/Ds_pp_MC_Prompt_Phi_pt${ptbins}.root" ,"/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_offMC/DsMinTree_pp_offMC_Prompt_Phi_pt2.root", 0 ,0)'
#	echo "/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/offpp_MC/Ds_pp_MC_Prompt_Phi_pt${ptbins[i]}.root"
	./MakeDsMinTree.exe "/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/offpp_MC/Ds_pp_MC_Prompt_Phi_pt${ptbins[i]}.root" "/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_offMC/DsMinTree_pp_offMC_Prompt_Phi_pt${ptbins[i]}.root" 0 0
	./MakeDsMinTree.exe "/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/offpp_MC/Ds_pp_MC_NonPrompt_Phi_pt${ptbins[i]}.root" "/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_offMC/DsMinTree_pp_offMC_NonPrompt_Phi_pt${ptbins[i]}.root" 0 0
	./MakeDsMinTree.exe "/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/offpp_MC/Ds_pp_MC_Prompt_f0_pt${ptbins[i]}.root" "/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_offMC/DsMinTree_pp_offMC_Prompt_f0_pt${ptbins[i]}.root" 0 0
	./MakeDsMinTree.exe "/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/offpp_MC/Ds_pp_MC_NonPrompt_f0_pt${ptbins[i]}.root" "/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_offMC/DsMinTree_pp_offMC_NonPrompt_f0_pt${ptbins[i]}.root" 0 0

done



#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/offpp_MC/DsNtuple_pp_Prompt_Phi_pt0.root" ,"output/DsMinTree_pp_offMC_Ds_Prompt_phikkpi_pt0.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/offpp_MC/Ds_pp_MC_Prompt_Phi_pt2.root" ,"/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_offMC/DsMinTree_pp_offMC_Prompt_Phi_pt2.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/offpp_MC/DsNtuple_pp_Prompt_Phi_pt4.root" ,"output/DsMinTree_pp_offMC_Ds_Prompt_phikkpi_pt4.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/offpp_MC/DsNtuple_pp_Prompt_Phi_pt6.root" ,"output/DsMinTree_pp_offMC_Ds_Prompt_phikkpi_pt6.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/offpp_MC/DsNtuple_pp_Prompt_Phi_pt10.root" ,"output/DsMinTree_pp_offMC_Ds_Prompt_phikkpi_pt10.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/offpp_MC/DsNtuple_pp_Prompt_Phi_pt19.root" ,"output/DsMinTree_pp_offMC_Ds_Prompt_phikkpi_pt19.root", 0 ,0)'


#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_Prompt_phikkpi_pt0.root" ,"output/DsMinTree_pp_MC_Ds_Prompt_phikkpi_pt0.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_Prompt_phikkpi_pt1p8.root" ,"output/DsMinTree_pp_MC_Ds_Prompt_phikkpi_pt1p8.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_Prompt_phikkpi_pt3p8.root" ,"output/DsMinTree_pp_MC_Ds_Prompt_phikkpi_pt3p8.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_Prompt_phikkpi_pt5p7.root" ,"output/DsMinTree_pp_MC_Ds_Prompt_phikkpi_pt5p7.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_Prompt_phikkpi_pt9p5.root" ,"output/DsMinTree_pp_MC_Ds_Prompt_phikkpi_pt9p5.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_Prompt_phikkpi_pt19.root" ,"output/DsMinTree_pp_MC_Ds_Prompt_phikkpi_pt19.root", 0 ,0)'

#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_NonPrompt_phikkpi_pt0.root" ,"output/DsMinTree_pp_MC_Ds_NonPrompt_phikkpi_pt0.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_NonPrompt_phikkpi_pt1p8.root" ,"output/DsMinTree_pp_MC_Ds_NonPrompt_phikkpi_pt1p8.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_NonPrompt_phikkpi_pt3p8.root" ,"output/DsMinTree_pp_MC_Ds_NonPrompt_phikkpi_pt3p8.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_NonPrompt_phikkpi_pt5p7.root" ,"output/DsMinTree_pp_MC_Ds_NonPrompt_phikkpi_pt5p7.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_NonPrompt_phikkpi_pt9p5.root" ,"output/DsMinTree_pp_MC_Ds_NonPrompt_phikkpi_pt9p5.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_NonPrompt_phikkpi_pt19.root" ,"output/DsMinTree_pp_MC_Ds_NonPrompt_phikkpi_pt19.root", 0 ,0)'
#
#
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_Prompt_f0980kkpi_pt0.root" ,"output/DsMinTree_pp_MC_Ds_Prompt_f0980kkpi_pt0.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_Prompt_f0980kkpi_pt1p8.root" ,"output/DsMinTree_pp_MC_Ds_Prompt_f0980kkpi_pt1p8.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_Prompt_f0980kkpi_pt3p8.root" ,"output/DsMinTree_pp_MC_Ds_Prompt_f0980kkpi_pt3p8.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_Prompt_f0980kkpi_pt5p7.root" ,"output/DsMinTree_pp_MC_Ds_Prompt_f0980kkpi_pt5p7.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_Prompt_f0980kkpi_pt9p5.root" ,"output/DsMinTree_pp_MC_Ds_Prompt_f0980kkpi_pt9p5.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_Prompt_f0980kkpi_pt19.root" ,"output/DsMinTree_pp_MC_Ds_Prompt_f0980kkpi_pt19.root", 0 ,0)'
#
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_NonPrompt_f0980kkpi_pt0.root" ,"output/DsMinTree_pp_MC_Ds_NonPrompt_f0980kkpi_pt0.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_NonPrompt_f0980kkpi_pt1p8.root" ,"output/DsMinTree_pp_MC_Ds_NonPrompt_f0980kkpi_pt1p8.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_NonPrompt_f0980kkpi_pt3p8.root" ,"output/DsMinTree_pp_MC_Ds_NonPrompt_f0980kkpi_pt3p8.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_NonPrompt_f0980kkpi_pt5p7.root" ,"output/DsMinTree_pp_MC_Ds_NonPrompt_f0980kkpi_pt5p7.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_NonPrompt_f0980kkpi_pt9p5.root" ,"output/DsMinTree_pp_MC_Ds_NonPrompt_f0980kkpi_pt9p5.root", 0 ,0)'
#root -b -l -q 'MakeDsMinTree.C++ ("/scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/pp_MC/Dntuple_pp_MC_NonPrompt_f0980kkpi_pt19.root" ,"output/DsMinTree_pp_MC_Ds_NonPrompt_f0980kkpi_pt19.root", 0 ,0)'
#
