#!/bin/bash
# run from within repo 

EXPERIMENT=$1

codeHomeDir=$2
PAIRED_OR_SINGLE_END_FRAGMENTS=$3
genomeDirSTAR=$4

STARFLAGS=${@:5} # pass all arguments after the first four

CURRENTEXPNUMBER=1
CURRENTSAMPLENUMBER=1
for fileName in $EXPERIMENT/raw/* ; do
    echo "Working on sample $CURRENTSAMPLENUMBER of experiment $CURRENTEXPNUMBER"
    SAMPLEID=`basename "$fileName"`

    fullCmd="$codeHomeDir/rajlabseqtools/Utilities/stepTwoStar/runStar.sh $EXPERIMENT $SAMPLEID $PAIRED_OR_SINGLE_END_FRAGMENTS $genomeDirSTAR $STARFLAGS"
    echo "$fullCmd"
    bsub -M "32000" -J "$EXPERIMENT.STAR.$CURRENTSAMPLENUMBER" -o "$EXPERIMENT/analyzed/$SAMPLEID/log/$(date +%Y-%m-%d_%H-%M).star.bsub.stdout" -e "$EXPERIMENT/analyzed/$SAMPLEID/log/$(date +%Y-%m-%d_%H-%M).star.bsub.stderr" -n 4 "$fullCmd"
    echo ""

    CURRENTSAMPLENUMBER=$((CURRENTSAMPLENUMBER+1))
done
echo "Submitted all samples to STAR"