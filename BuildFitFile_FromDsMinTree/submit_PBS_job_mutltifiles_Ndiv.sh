#!/bin/sh

if [ "$#" -ne 7 ]; then
	echo "Wrong number of parameters. 7 expected, passed $#"
	exit
fi

INFILELIST=$1
isPbPb=$2  # need to be set in loop.C
FilePerJob=$3
OUTPUTDIR=$4
JobIndex=$5
isReal=$6
Ndiv=$7
TotalFiles=$(wc -l < "$INFILELIST")
Njobs=$((TotalFiles/FilePerJob+1))

DIR="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/BuildFitFile_FromDsMinTree"

echo "INFILELIST: $INFILELIST"
echo "Total number of file: $TotalFiles, FilePerJob: $FilePerJob"
echo "$Njobs jobs will be submitted"

cd ${DIR}



if [ ! -d $OUTPUTDIR ]; then
    mkdir -p $OUTPUTDIR
fi

#mkdir -p PBS_jobs/PBS_log/${JobIndex}
#mkdir -p PBS_jobs/PBS_shAll/${JobIndex}

mkdir -p PBS_jobs/PBS_log/
mkdir -p PBS_jobs/PBS_shAll/

for ((count=1; count <= $Njobs; count++))
do
	echo "Job $count"

	#decide start file and end file
	STARTFILE=$(((count-1)*FilePerJob))
	ENDFILE=$((count*FilePerJob))
	if [ $ENDFILE -gt $TotalFiles ]; then
		ENDFILE=$TotalFiles
	fi
	echo "Files from $STARTFILE to $ENDFILE"

	for ((idiv=0; idiv < $Ndiv; idiv++))
	do
# make the PBS file
cat > PBS_jobs/PBS_shAll/PBS_${STARTFILE}to${ENDFILE}_${JobIndex}_${idiv}.sh <<EOF
#!/bin/sh

#PBS -l nodes=1,walltime=04:00:00,mem=2400MB
#PBS -N Dnt${STARTFILE}to${ENDFILE}_${JobIndex}_${idiv}
#PBS -o ${DIR}/PBS_jobs/PBS_log/PBS_${STARTFILE}to${ENDFILE}_${JobIndex}_${idiv}.log
#PBS -e ${DIR}/PBS_jobs/PBS_log/PBS_${STARTFILE}to${ENDFILE}_${JobIndex}_${idiv}.err
#PBS -r n
#PBS -V
#PBS -q standby 

cd ${DIR}

./exec_condorfilelist_Ndiv.sh $INFILELIST $isPbPb $STARTFILE $ENDFILE $OUTPUTDIR $isReal $Ndiv $idiv

EOF

	chmod 744 PBS_jobs/PBS_shAll/PBS_${STARTFILE}to${ENDFILE}_${JobIndex}_${idiv}.sh

	qsub PBS_jobs/PBS_shAll/PBS_${STARTFILE}to${ENDFILE}_${JobIndex}_${idiv}.sh
	echo "PBS_jobs/PBS_shAll/PBS_${STARTFILE}to${ENDFILE}_${JobIndex}_${idiv}.sh submitted"

	done #done idiv

done #done count



##/home/peng43/work/Project/Ds_PbPb/CMSSW/DsFinder/working_DsFinder/CMSSW_7_5_8_patch3/src/Bfinder/Bfinder/Dntuple/PBS_jobs/exec_condorfilelist.sh $INFILELIST $isPbPb $STARTFILE $ENDFILE $OUTPUTDIR $isReal
