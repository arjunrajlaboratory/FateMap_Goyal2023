#!
# run from within "repo" directory

EXPERIMENT=$1
SAMPLEID=$2
PAIRED_OR_SINGLE_END_FRAGMENTS=$3
genomeDirSTAR=$4
STARFLAGS=${@:5} # pass all arguments after the first two

toolNAME=star

if [ ! -d $EXPERIMENT/analyzed ]; then
    mkdir $EXPERIMENT/analyzed
fi

if [ ! -d $EXPERIMENT/analyzed/$SAMPLEID ]; then
    mkdir $EXPERIMENT/analyzed/$SAMPLEID
fi

if [ ! -d $EXPERIMENT/analyzed/$SAMPLEID/log ]; then
    mkdir $EXPERIMENT/analyzed/$SAMPLEID/log
fi


destinationDir="$EXPERIMENT/analyzed/$SAMPLEID/star/"
if [ ! -d $destinationDir ]; then
    mkdir "$destinationDir"
fi

inputFileR1="$EXPERIMENT/raw/$SAMPLEID/${SAMPLEID}_R1.fastq.gz"
inputFileR2="$EXPERIMENT/raw/$SAMPLEID/${SAMPLEID}_R2.fastq.gz"

numCPU=4
if [ $PAIRED_OR_SINGLE_END_FRAGMENTS = "single" ]; then
	cmdToRun="STAR \
		--genomeDir $genomeDirSTAR \
		--readFilesIn $inputFileR1 \
		--genomeLoad LoadAndRemove \
		--outFileNamePrefix $destinationDir/$SAMPLEID. \
		--runThreadN $numCPU \
		$STARFLAGS"
fi

if [ $PAIRED_OR_SINGLE_END_FRAGMENTS = "paired" ]; then
	cmdToRun="STAR \
		--genomeDir $genomeDirSTAR \
		--readFilesIn $inputFileR1 $inputFileR2 \
		--genomeLoad LoadAndRemove \
		--outFileNamePrefix $destinationDir/$SAMPLEID. \
		--runThreadN $numCPU \
		$STARFLAGS"
fi

JOURNAL=$EXPERIMENT/analyzed/$SAMPLEID/log/$(date +%Y-%m-%d_%H-%M).$toolNAME.log
echo "Starting..." >> $JOURNAL
date >> $JOURNAL
echo "$cmdToRun" >> $JOURNAL
eval "$cmdToRun"
date >> $JOURNAL
echo "Done" >> $JOURNAL
