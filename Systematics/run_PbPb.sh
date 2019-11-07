#!/bin/sh

# root -b -q 'CutScanSys.C++(3)' > ./log/custScanSys_PbPb.log
# root -b -q 'CutScanSys_FixShape.C++(3)' > ./log/custScanSys_PbPb_FixShape.log
root -b -q 'BtoDsSys.C++(3)' > ./log/BtoDsSys_PbPb.log
root -b -q 'MCShape_f0Shape_Sys.C++(3)' > ./log/MCShape_PbPb.log
root -b -q 'SystematicSum.C++(3)' > ./log/SysSum_PbPb.log


