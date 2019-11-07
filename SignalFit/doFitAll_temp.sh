#!/bin/sh

CMSenvNew

echo "do fit pp"
# ./doFit_pp.sh
./doFit_pp_temp.sh
echo "fit pp done"
echo "do fit PbPB"
# ./doFit_PbPb3.sh
./doFit_PbPb3_temp.sh
echo "fit PbPb done"
