### for PbPb_GJ


gfal-copy -f -t 50000000 -T 25000000 file://///scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/GJ/   gsiftp://cms-gridftp.rcac.purdue.edu/store/group/hi/chengchi/Dntuple/DsMinTree/PbPb_Data/GJ/ -f  2>&1 | tee DsMinTree_GJ.log & 

gfal-copy -f -t 50000000 -T 25000000 file://///scratch/halstead/p/peng43/Ds_phikkpi/DsMinTree/PbPb_Data/TrkOnly/   gsiftp://cms-gridftp.rcac.purdue.edu/store/group/hi/chengchi/Dntuple/DsMinTree/PbPb_Data/TrkOnly/ -f  2>&1 | tee DsMinTree_TrkOnly.log & 

#### for PbPb_Trkonly


