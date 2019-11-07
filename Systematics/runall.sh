#!/bin/sh

./run_pp.sh

./run_PbPb.sh

root -b -q SystematicSumRaa.C++> ./log/SysRaaSum.log
