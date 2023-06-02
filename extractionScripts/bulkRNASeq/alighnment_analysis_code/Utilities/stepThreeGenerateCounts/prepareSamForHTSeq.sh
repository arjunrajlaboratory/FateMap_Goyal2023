#!
# run from within "repo" directory

EXPERIMENT=$1
SAMPLEID=$2

ALIGNMENT_TOOL_NAME=rum

if [[ ! -z "$3" ]]; then
    ALIGNMENT_TOOL_NAME="$3"
fi

commandName=prepareSamForHTSeq

if [ ! -d $EXPERIMENT/analyzed/$SAMPLEID/log ]; then
    mkdir $EXPERIMENT/analyzed/$SAMPLEID/log
fi
JOURNAL=$EXPERIMENT/analyzed/$SAMPLEID/log/$(date +%Y-%m-%d_%H-%M).$commandName.log


if [ ! -d $EXPERIMENT/analyzed/$SAMPLEID/htseq ]; then
    mkdir $EXPERIMENT/analyzed/$SAMPLEID/htseq
fi

readNameSortCmd="samtools sort -n \
    $EXPERIMENT/analyzed/$SAMPLEID/$ALIGNMENT_TOOL_NAME/$SAMPLEID.sorted.mapped.unique.bam \
    $EXPERIMENT/analyzed/$SAMPLEID/htseq/$SAMPLEID.nameSorted"

fixMateCmd="samtools fixmate \
    $EXPERIMENT/analyzed/$SAMPLEID/htseq/$SAMPLEID.nameSorted.bam \
    $EXPERIMENT/analyzed/$SAMPLEID/htseq/$SAMPLEID.nameSorted.mateFixed.bam"

toSamCmd="samtools view -h \
    $EXPERIMENT/analyzed/$SAMPLEID/htseq/$SAMPLEID.nameSorted.mateFixed.bam \
    > $EXPERIMENT/analyzed/$SAMPLEID/htseq/$SAMPLEID.nameSorted.mateFixed.sam"

echo "Starting..." >> $JOURNAL
date >> $JOURNAL
echo "$readNameSortCmd" >> $JOURNAL
eval "$readNameSortCmd"
echo "$fixMateCmd" >> $JOURNAL
eval "$fixMateCmd"
echo "$toSamCmd" >> $JOURNAL
eval "$toSamCmd"
date >> $JOURNAL
echo "Done" >> $JOURNAL
