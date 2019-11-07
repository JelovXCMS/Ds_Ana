#!/bin/sh

INPUTFILE=$1

gfal-copy file:////scratch/halstead/p/peng43/Ds_phikkpi/${INPUTFILE} gsiftp://cms-gridftp.rcac.purdue.edu/store/user/chengchi/DsNtuple/${INPUTFILE}
