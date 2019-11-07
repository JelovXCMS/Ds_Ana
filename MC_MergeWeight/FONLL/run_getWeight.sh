#!/bin/sh

root -b -q 'getFONLLweight.C++(0,0,0)' > log/pp_Prompt_phi.log
#root -b -q 'getFONLLweight.C++(0,0,1)' > log/pp_Prompt_f0.log
root -b -q 'getFONLLweight.C++(0,1,0)' > log/pp_NonPrompt_phi.log
#root -b -q 'getFONLLweight.C++(0,1,1)' > log/pp_NonPrompt_f0.log

echo " pp done"

root -b -q 'getFONLLweight.C++(1,0,0)' > log/PbPb_Prompt_phi.log
#root -b -q 'getFONLLweight.C++(1,0,1)' > log/PbPb_Prompt_f0.log
root -b -q 'getFONLLweight.C++(1,1,0)' > log/PbPb_NonPrompt_phi.log
#root -b -q 'getFONLLweight.C++(1,1,1)' > log/PbPb_NonPrompt_f0.log

echo " PbPb done"
