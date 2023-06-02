#!
# run from within "repo" directory

EXPERIMENT=$1
SAMPLEID=$2
STARFLAGS=${@:3} # pass all arguments after the first two

toolNAME=bowtie2

if [ ! -d $EXPERIMENT/analyzed ]; then
    mkdir $EXPERIMENT/analyzed
fi

if [ ! -d $EXPERIMENT/analyzed/$SAMPLEID ]; then
    mkdir $EXPERIMENT/analyzed/$SAMPLEID
fi

if [ ! -d $EXPERIMENT/analyzed/$SAMPLEID/log ]; then
    mkdir $EXPERIMENT/analyzed/$SAMPLEID/log
fi


destinationDir="$EXPERIMENT/analyzed/$SAMPLEID/bowtie2"
if [ ! -d $destinationDir ]; then
    mkdir "$destinationDir"
fi

inputFile="$EXPERIMENT/raw/$SAMPLEID/$SAMPLEID.fastq"
if [ ! -e $inputFile ]; then
    gunzip -c "$inputFile.gz" > "$inputFile"
fi

# bowtie2 location
bowtie2NAME="/Users/sydneyshaffer/ATACseqTools/bowtie2-2.2.6/bowtie2"

# this is the index downloaded from the bowtie2 website
genomeDir="/Users/sydneyshaffer/ATACseqTools/hg19/hg19"

numCPU=3
cmdToRun="bowtie2 -x $genomeDir -U $inputFile -S $destinationDir/$SAMPLEID.sam -p $numCPU"

JOURNAL=$EXPERIMENT/analyzed/$SAMPLEID/log/$(date +%Y-%m-%d_%H-%M).$toolNAME.log
echo "Starting..." >> $JOURNAL
date >> $JOURNAL
echo "$cmdToRun" >> $JOURNAL
eval "$cmdToRun"
date >> $JOURNAL
echo "Done" >> $JOURNAL
