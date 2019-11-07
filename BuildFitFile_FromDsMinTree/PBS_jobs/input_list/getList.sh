
ls /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/GJ/Ds_PbPb_Data_HIMB2_GJ/*.root >  Ds_PbPb_Data_HIMB2_GJ.lis
ls /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/GJ/Ds_PbPb_Data_HIMB3_GJ/*.root >  Ds_PbPb_Data_HIMB3_GJ.lis
ls /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/GJ/Ds_PbPb_Data_HIMB4_GJ/*.root >  Ds_PbPb_Data_HIMB4_GJ.lis


#ls /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/TrkOnly/Ds_PbPb_Data_HIMB1_TrkOnly/*.root >  Ds_PbPb_Data_HIMB_TrkOnly.lis

for i in {1..11}
	do
	ls /scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/TrkOnly/Ds_PbPb_Data_HIMB${i}_TrkOnly/*.root >  Ds_PbPb_Data_HIMB${i}_TrkOnly.lis
done
