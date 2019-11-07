#!/bin/sh

if [ "$#" -ne 6 ]; then
	echo "Wrong number of parameters. 5 expected, passed $#"
	exit
fi

INFILELIST=$1
isPbPb=$2  # need to be set in loop.C
FilePerJob=$3
OUTPUTDIR=$4
JobIndex=$5
isReal=$6
TotalFiles=$(wc -l < "$INFILELIST")
Njobs=$((TotalFiles/FilePerJob+1))

DIR="/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/MC_CrossSection"

echo "INFILELIST: $INFILELIST"
echo "Total number of file: $TotalFiles, FilePerJob: $FilePerJob"
echo "$Njobs jobs will be submitted"

cd ${DIR}


#rm MakeDsMinTree.exe
#g++ MakeDsMinTree.C $(root-config --cflags --libs) -Wall -O2 -o MakeDsMinTree.exe

## put other needed exe here


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

# make the PBS file
cat > PBS_jobs/PBS_shAll/PBS_${STARTFILE}to${ENDFILE}_${JobIndex}.sh <<EOF
#!/bin/sh

#PBS -l nodes=1,walltime=04:00:00,mem=3600MB
#PBS -N Dnt${STARTFILE}to${ENDFILE}_${JobIndex}
#PBS -o ${DIR}/PBS_jobs/PBS_log/PBS_${STARTFILE}to${ENDFILE}_${JobIndex}.log
#PBS -e ${DIR}/PBS_jobs/PBS_log/PBS_${STARTFILE}to${ENDFILE}_${JobIndex}.err
#PBS -r n
#PBS -V
#PBS -q standby 

cd ${DIR}

./exec_condorfilelist.sh $INFILELIST $isPbPb $STARTFILE $ENDFILE $OUTPUTDIR $isReal

EOF

	chmod 744 PBS_jobs/PBS_shAll/PBS_${STARTFILE}to${ENDFILE}_${JobIndex}.sh

	qsub PBS_jobs/PBS_shAll/PBS_${STARTFILE}to${ENDFILE}_${JobIndex}.sh
	echo "PBS_jobs/PBS_shAll/PBS_${STARTFILE}to${ENDFILE}_${JobIndex}.sh submitted"

done



##/home/peng43/work/Project/Ds_PbPb/CMSSW/DsFinder/working_DsFinder/CMSSW_7_5_8_patch3/src/Bfinder/Bfinder/Dntuple/PBS_jobs/exec_condorfilelist.sh $INFILELIST $isPbPb $STARTFILE $ENDFILE $OUTPUTDIR $isReal
