#!/bin/sh

CMSenvNew

echo "do fit pp"
source ./doFit_pp.sh
echo "fit pp done"
echo "do fit PbPB"
source ./doFit_PbPb3.sh
echo "fit PbPb done"
