
# int FitCutScan(int isPbPb=0, double DptLow=2, double DptHigh=6, TString var_scan="Ddls", double scan_val=2.5, TString scanType="larger")

isPbPb=0
DptLow=2
DptHigh=6
var_scan="Dchi2cl"
# scan_val=2.5

# scan_val=(1.5 1.75 2.0 2.25 2.5 2.75 3.0 3.25 3.5 3.75 4.0)
# scan_val=(3.0 3.5 4.0 4.5 5.0)
# scan_val=(2.75 3.0 3.25 3.5 3.75 4.0 4.25 4.5 4.75 5.0 5.25 5.5)
scan_val=(0.02 0.07 0.12 0.17 0.22 0.27 0.32)
# scan_val=(0.05 0.10 0.15 0.20 0.25 0.3 0.35 0.4)
# scan_val=(0.45 0.50)
# scan_val=(0.05 0.1 0.15 0.20 0.25 0.3 0.35 0.4)
# scan_val=(0.175 0.225 0.275 0.325 0.375 0.425)
# scan_val=(0.2 0.25 0.3 0.35 )
# scan_val=(0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4)
# scan_val=(0.15 0.45)
# scan_val=(2.0 3.5)
FixShape=1
FixShapeVal=0.17


jobscount=0

rm FitCutScan.exe
# g++ FitCutScan.C $(root-config --cflags --libs) -L $ROOTSYS/lib -lRooFit -lHtml -lMinuit -lRooFitCore -Wall -O2 -o FitCutScan.exe
g++ FitCutScan.C $(root-config --cflags --libs) -L $ROOTSYS/lib -lRooFit -lHtml -lMinuit -lRooFitCore -lRooStats   -Wall -O2 -o FitCutScan.exe
# g++ FitCutScan.C $(root-config --cflags --libs)  $(root-config --ldflags --glibs ) -L $ROOTSYS/lib -lRooFit -lHtml -lMinuit -lRooFitCore -Wall -O2 -o FitCutScan.exe

for i in ${scan_val[@]}
do
	echo ${var_scan} " = " ${i}

	./FitCutScan.exe  ${isPbPb}  ${DptLow}  ${DptHigh}  ${var_scan}  ${i}  ${FixShape} ${FixShapeVal} >log/Fit_isPbPb${isPbPb}_Dpt${DptLow}to${DptHigh}_${Ddls}${i}_FixShape${FixShape}.log &

	jobscount=$((${jobscount}+1))
	echo "jobscount = "${jobscount}
	if [ $((${jobscount} % 4)) -eq 0 ]
	then
		echo "wait"
		wait
	fi


done

wait
# root -b -q 'FitCutScan(${isPbPb},)'
