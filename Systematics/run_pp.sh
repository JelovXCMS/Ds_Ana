#!/bin/sh

# root -b -q 'CutScanSys.C++(0)' > ./log/custScanSys_pp.log

# root -b -q 'CutScanSys_FixShape.C++(0)' > ./log/custScanSys_pp_FixShape.log // replaced by ARC_HW/fitCutscan

root -b -q 'BtoDsSys.C++(0)' > ./log/BtoDsSys_pp.log
root -b -q 'MCShape_f0Shape_Sys.C++(0)' > ./log/MCShape_pp.log
root -b -q 'SystematicSum.C++(0)' > ./log/SysSum_pp.log


