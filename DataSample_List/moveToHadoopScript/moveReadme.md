example :



screen -S gfal-copy


"Ctrl-a" followed by "d"
That will disconnect you from the screen session, but
gfal-copy will not know about that.

You can now close your SSH connection to hep.rcac
When you login later, resume your screen-session by doing:
'screen -dr mylongsessionname'


gfal-copy -r -f -t 50000000 -T 25000000 file://///scratch/halstead/p/peng43/Ds_phikkpi/Dntuple gsiftp://cms-gridftp.rcac.purdue.edu/store/group/hi/chengchi/Ds_phikkpi_0510bck/Dntuple  -f 2>&1 | tee DsNtuple_toMNT.log

gfal-copy -r -f -t 50000000 -T 25000000 file://///scratch/halstead/p/peng43/Ds_phikkpi/Dntuple/PbPb_Data/GJ/HIMB2_GJpart1/   gsiftp://cms-gridftp.rcac.purdue.edu/store/group/hi/chengchi/Dntuple/Dsntuple/PbPb_Data/GJ/HIMB2_GJpart1/ -f  2>&1 | tee DsNtuple_toMNT.log


