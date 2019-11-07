#!/bin/sh

# ./run_pp.sh

# ./run_PbPb.sh

root -b -q 'SystematicSum.C++(0)' > ./log/SysSum_pp.log
root -b -q 'SystematicSum.C++(3)' > ./log/SysSum_PbPb.log

root -b -q SystematicSumRaa.C++> ./log/SysRaaSum.log
