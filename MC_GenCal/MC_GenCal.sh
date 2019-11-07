#!/bin/sh

#	echo example usage : source  MC_GenCal.sh Ds_Gen_Prompt_phikkpi_pt1p8 /mnt/hadoop/store/user/chengchi/Ds_ppMC_180223/Ds_Prompt_phikkpi_pp/Ds_Gen_Prompt_phikkpi_pt1p8_Pythia8_180223/180303_145314/0000/log

 Type=$1
 Dir=$2

 echo $Type $Dir
	
	mkdir -p Summary
	mkdir -p $Type
	cd $Type

	tar zxvf ${Dir}/0000/log/cmsRun_3.log.tar.gz
	tar zxvf ${Dir}/0000/log/cmsRun_5.log.tar.gz
	tar zxvf ${Dir}/0000/log/cmsRun_7.log.tar.gz
	tar zxvf ${Dir}/0000/log/cmsRun_14.log.tar.gz
	tar zxvf ${Dir}/0000/log/cmsRun_18.log.tar.gz

	du -BM ${Dir} > Summary_${Type}.txt

	grep "Filter efficiency"  cmsRun-stdout-3.log  >>Summary_${Type}.txt
	grep "cross section ="  cmsRun-stdout-3.log  >>Summary_${Type}.txt
	grep "Avg event"  cmsRun-stdout-3.log  >>Summary_${Type}.txt
	grep "Total job"  cmsRun-stdout-3.log  >>Summary_${Type}.txt

	grep "Filter efficiency"  cmsRun-stdout-5.log  >>Summary_${Type}.txt
	grep "cross section ="  cmsRun-stdout-5.log  >>Summary_${Type}.txt
	grep "Avg event"  cmsRun-stdout-5.log  >>Summary_${Type}.txt
	grep "Total job"  cmsRun-stdout-5.log  >>Summary_${Type}.txt


	grep "Filter efficiency"  cmsRun-stdout-7.log  >>Summary_${Type}.txt
	grep "cross section ="  cmsRun-stdout-7.log  >>Summary_${Type}.txt
	grep "Avg event"  cmsRun-stdout-7.log  >>Summary_${Type}.txt
	grep "Total job"  cmsRun-stdout-7.log  >>Summary_${Type}.txt

	grep "Filter efficiency"  cmsRun-stdout-14.log  >>Summary_${Type}.txt
	grep "cross section ="  cmsRun-stdout-14.log  >>Summary_${Type}.txt
	grep "Avg event"  cmsRun-stdout-14.log  >>Summary_${Type}.txt
	grep "Total job"  cmsRun-stdout-14.log  >>Summary_${Type}.txt

	grep "Filter efficiency"  cmsRun-stdout-18.log  >>Summary_${Type}.txt
	grep "cross section ="  cmsRun-stdout-18.log  >>Summary_${Type}.txt
	grep "Avg event"  cmsRun-stdout-18.log  >>Summary_${Type}.txt
	grep "Total job"  cmsRun-stdout-18.log  >>Summary_${Type}.txt

	mv Summary_${Type}.txt ../Summary

	rm *.log
  rm *.xml
	cd ..
	rm -r $Type

