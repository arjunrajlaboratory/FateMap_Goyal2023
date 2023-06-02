#!
# run from within "repo" directory

EXPERIMENT=$1
SAMPLEID=$2

ALIGNMENT_TOOL_NAME=bowtie2

COMMANDNAME=HOMER

if [ ! -d $EXPERIMENT/analyzed/$SAMPLEID/HOMER ]; then
    mkdir $EXPERIMENT/analyzed/$SAMPLEID/HOMER
fi

if [ ! -d $EXPERIMENT/analyzed/$SAMPLEID/log ]; then
    mkdir $EXPERIMENT/analyzed/$SAMPLEID/log
fi

JOURNAL=$EXPERIMENT/analyzed/$SAMPLEID/log/$(date +%Y-%m-%d_%H-%M).$COMMANDNAME.log


# runs HOMER
COMMAND1="makeTagDirectory $EXPERIMENT/analyzed/$SAMPLEID/HOMER \
	$EXPERIMENT/analyzed/$SAMPLEID/$SAMPLEID.sorted.mapped.rmDups.rmMandY.bam"

COMMAND2="makeUCSCfile $EXPERIMENT/analyzed/$SAMPLEID/HOMER -o auto"

# -o auto means this will be written to the tag directory (here just HOMER/)
COMMAND3="findPeaks $EXPERIMENT/analyzed/$SAMPLEID/HOMER -style dnase -o auto \
	> $EXPERIMENT/analyzed/$SAMPLEID/HOMER/peaks.txt"
	
COMMAND4="pos2bed.pl $EXPERIMENT/analyzed/$SAMPLEID/HOMER/peaks.txt \
	> $EXPERIMENT/analyzed/$SAMPLEID/HOMER/peaks.bed"

    
echo "Starting..." >> $JOURNAL

date >> $JOURNAL
#echo $COMMAND1 >> $JOURNAL
#eval $COMMAND1
 
#echo $COMMAND2 >> $JOURNAL
#eval $COMMAND2
 
echo $COMMAND3 >> $JOURNAL
eval $COMMAND3

echo $COMMAND4 >> $JOURNAL
eval $COMMAND4