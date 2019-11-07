#!/bin/sh

gfal-copy -f -t 50000000 -T 25000000 file:////mnt/hadoop/store/user/jisun/Dsfinder_phikkpi_180305/PbPb_data/MB_TrkOnly   gsiftp://cms-gridftp.rcac.purdue.edu/store/group/hi/chengchi/Dsfinder_phikkpi_180305/PbPb_data/MB_TrkOnly 2>&1 | tee copy_Dsfinder_PbPb_MB_Trkonly_fromJian.log
