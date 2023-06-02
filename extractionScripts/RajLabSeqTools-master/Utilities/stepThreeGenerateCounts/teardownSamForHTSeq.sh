#!
# run from within "repo" directory

EXPERIMENT=$1
SAMPLEID=$2

commandName=teardownSamForHTSeq

if [ ! -d $EXPERIMENT/analyzed/$SAMPLEID/log ]; then
    mkdir $EXPERIMENT/analyzed/$SAMPLEID/log
fi
JOURNAL=$EXPERIMENT/analyzed/$SAMPLEID/log/$(date +%Y-%m-%d_%H-%M).$commandName.log

nameSortedBamFile="$EXPERIMENT/analyzed/$SAMPLEID/htseq/$SAMPLEID.nameSorted.bam"
mateFixedBamFile="$EXPERIMENT/analyzed/$SAMPLEID/htseq/$SAMPLEID.nameSorted.mateFixed.bam"
mateFixedSamFile="$EXPERIMENT/analyzed/$SAMPLEID/htseq/$SAMPLEID.nameSorted.mateFixed.sam"

filesToLookFor=($nameSortedBamFile $mateFixedBamFile $mateFixedSamFile)

echo "Starting..." >> $JOURNAL
date >> $JOURNAL

for fileName in "${filesToLookFor[@]}" ; do
	if [ -e $fileName ]; then
		cmdToRun="unlink $fileName"
		echo "$cmdToRun" >> $JOURNAL
		eval "$cmdToRun"
	fi
done

date >> $JOURNAL
echo "Done" >> $JOURNAL
