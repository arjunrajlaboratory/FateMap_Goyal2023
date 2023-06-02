#!
# run from within "repo" directory

EXPERIMENT=$1
SAMPLEID=$2

ALIGNMENT_TOOL_NAME=bowtie2

COMMANDNAME=macs2

if [ ! -d $EXPERIMENT/analyzed/$SAMPLEID/MACS2 ]; then
    mkdir $EXPERIMENT/analyzed/$SAMPLEID/MACS2
fi

JOURNAL=$EXPERIMENT/analyzed/$SAMPLEID/log/$(date +%Y-%m-%d_%H-%M).$COMMANDNAME.log

MACS2dir="/Users/sydneyshaffer/ATACseqTools/MACS-master/MACS2"


# runs MACS2
COMMAND="macs2 callpeak -t \
	$EXPERIMENT/analyzed/$SAMPLEID/$ALIGNMENT_TOOL_NAME/$SAMPLEID.sorted.mapped.rmDups.rmMandY.bam \
	-n $SAMPLEID --outdir $EXPERIMENT/analyzed/$SAMPLEID/MACS2/ -f BAM \
	-g hs --nomodel -B --SPMR --nolambda --keep-dup all --call-summits"
    
echo "Starting..." >> $JOURNAL

date >> $JOURNAL
echo $COMMAND >> $JOURNAL
eval $COMMAND