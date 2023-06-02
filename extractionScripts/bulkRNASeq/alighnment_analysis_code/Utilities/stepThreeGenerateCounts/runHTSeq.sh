 #!
# run from within "repo" directory

EXPERIMENT=$1
SAMPLEID=$2
gtfFile=$3

commandNAME=runHTSeq

if [ ! -d $EXPERIMENT/analyzed/$SAMPLEID/log ]; then
    mkdir $EXPERIMENT/analyzed/$SAMPLEID/log
fi

JOURNAL=$EXPERIMENT/analyzed/$SAMPLEID/log/$(date +%Y-%m-%d_%H-%M).$commandNAME.log


if [ ! -d $EXPERIMENT/analyzed/$SAMPLEID/htseq ]; then
    mkdir $EXPERIMENT/analyzed/$SAMPLEID/htseq
fi

inputFile="$EXPERIMENT/analyzed/$SAMPLEID/htseq/$SAMPLEID.nameSorted.mateFixed.sam"

countsOutFile="$EXPERIMENT/analyzed/$SAMPLEID/htseq/$SAMPLEID.htseq.stdout"
logOutFile="$EXPERIMENT/analyzed/$SAMPLEID/htseq/$SAMPLEID.htseq.stderr"

htseqCommand="python -m HTSeq.scripts.count \
	--idattr=gene_id \
	--mode=union \
	--stranded=no \
	--type=exon \
	--order=name \
	$inputFile \
	$gtfFile \
	1> $countsOutFile \
	2> $logOutFile"

echo "Starting..." >> $JOURNAL
date >> $JOURNAL
echo "$htseqCommand" >> $JOURNAL
eval "$htseqCommand"

date >> $JOURNAL
echo "Done" >> $JOURNAL
