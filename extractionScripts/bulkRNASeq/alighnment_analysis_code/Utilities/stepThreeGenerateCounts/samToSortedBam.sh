#!
# run from within "repo" directory

EXPERIMENT=$1
SAMPLEID=$2

ALIGNMENT_TOOL_NAME=bowtie2

if [[ ! -z "$3" ]]; then
    ALIGNMENT_TOOL_NAME="$3"
fi

COMMANDNAME=samtools

if [ ! -d $EXPERIMENT/analyzed/$SAMPLEID/log ]; then
    mkdir $EXPERIMENT/analyzed/$SAMPLEID/log
fi

JOURNAL=$EXPERIMENT/analyzed/$SAMPLEID/log/$(date +%Y-%m-%d_%H-%M).$COMMANDNAME.log

COMMAND="samtools view -bS \
    $EXPERIMENT/analyzed/$SAMPLEID/$ALIGNMENT_TOOL_NAME/$SAMPLEID.Aligned.out.sam \
    | samtools sort - \
    $EXPERIMENT/analyzed/$SAMPLEID/$ALIGNMENT_TOOL_NAME/$SAMPLEID.sorted"


INDEXCOMMAND="samtools index \
    $EXPERIMENT/analyzed/$SAMPLEID/$ALIGNMENT_TOOL_NAME/$SAMPLEID.sorted.bam \
    $EXPERIMENT/analyzed/$SAMPLEID/$ALIGNMENT_TOOL_NAME/$SAMPLEID.sorted.bai"

REMOVE_SAM_COMMAND="rm $EXPERIMENT/analyzed/$SAMPLEID/$ALIGNMENT_TOOL_NAME/$SAMPLEID.Aligned.out.sam"

echo "Starting..." >> $JOURNAL

date >> $JOURNAL
echo $COMMAND >> $JOURNAL
eval $COMMAND

date >> $JOURNAL
echo $INDEXCOMMAND >> $JOURNAL
eval $INDEXCOMMAND

date >> $JOURNAL
echo $REMOVE_SAM_COMMAND >> $JOURNAL
eval $REMOVE_SAM_COMMAND

date >> $JOURNAL
echo "Done" >> $JOURNAL
