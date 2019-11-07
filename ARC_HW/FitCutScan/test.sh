#!/bin/bash

ibin_Dpt=1

scan_val=(1.5 1.75 2.0 2.25 2.5 2.75 3.0 3.25 3.5 3.75 4.0)
scan_val1=(0.1 0.2)


echo "ibin_Dpt = " ibin_Dpt
if [ $((ibin_Dpt)) -eq 1]
then 
	scan_val=scan_val1
	echo 
