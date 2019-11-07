
# int FitCutScan(int isPbPb=0, double DptLow=2, double DptHigh=6, TString var_scan="Ddls", double scan_val=2.5, TString scanType="larger")

isPbPb=0
DptLow=2
DptHigh=6
var_scan="Ddls"
scan_val=2.5

scan_val=(1.5 1.75 2.0 2.25 2.5 2.75 3.0 3.25 3.5 3.75 4.0)
# scan_val=(3.5 4.0 4.5 5.0)
# scan_val=(3.25 3.5 3.75 4.0 4.25 4.5 4.75 5.0 5.25 5.5)
# scan_val=(2.75 3.0 3.25 3.5 3.75 4.0 4.25 4.5 4.75 5.0 5.25 5.5)
# scan_val=(2.5 2.75 3.0)
# scan_val=(2.0 3.5)
FixShape=1
FixShapeVal=3.25


jobscount=0

rm FitCutScan.exe
# g++ FitCutScan.C $(root-config --cflags --libs) -L $ROOTSYS/lib -lRooFit -lHtml -lMinuit -lRooFitCore -Wall -O2 -o FitCutScan.exe
g++ FitCutScan.C $(root-config --cflags --libs) -L $ROOTSYS/lib -lRooFit -lHtml -lMinuit -lRooFitCore -lRooStats   -Wall -O2 -o FitCutScan.exe
# g++ FitCutScan.C $(root-config --cflags --libs)  $(root-config --ldflags --glibs ) -L $ROOTSYS/lib -lRooFit -lHtml -lMinuit -lRooFitCore -Wall -O2 -o FitCutScan.exe

for i in ${scan_val[@]}
do
	echo ${var_scan} " = " ${i}

	./FitCutScan.exe  ${isPbPb}  ${DptLow}  ${DptHigh}  ${var_scan}  ${i}  ${FixShape} ${FixShapeVal}  >log/Fit_isPbPb${isPbPb}_Dpt${DptLow}to${DptHigh}_${Ddls}${i}_FixShape${FixShape}.log &

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
