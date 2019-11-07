#!/bin/sh

# root -b -q -l 'BuildFitFile_FromDsMinTreeLoop.C++("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_Data/DsMinTree_pp_Data_MBAll.root" , "./output/FitFile_pp", true, 0)'

# root -b -q -l 'BuildFitFile_FromDsMinTreeLoop_moreScan.C++("/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/DataSample_List/DsMinTree/pp_Data/DsMinTree_pp_Data_MBAll.root" , "./output/FitFile_pp_moreScan", true, 0)'
root -b -q -l 'BuildFitFile_FromDsMinTreeLoop_moreScan.C++("/scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/pp_Data/DsMinTree_pp_Data_MBAll.root" , "./output/FitFile_pp_moreScan", true, 0)'
