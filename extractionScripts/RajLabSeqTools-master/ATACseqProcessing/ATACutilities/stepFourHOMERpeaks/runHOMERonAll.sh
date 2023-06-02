#!/bin/bash



# run from within repo 

# runs HOMER on each file

EXPERIMENT=$1

codeHomeDir=$2


CURRENTEXPNUMBER=1
CURRENTSAMPLENUMBER=1
for fileName in $EXPERIMENT/analyzed/* ; do
    echo "Working on sample $CURRENTSAMPLENUMBER of experiment $CURRENTEXPNUMBER"
    echo "$sampleID"
    echo "$fileName"
    sampleID=`basename "$fileName"`

	# 
    fullCmd="$codeHomeDir/rajlabseqtools/ATACseqProcessing/ATACutilities/stepFourHOMERpeaks/runHOMER.sh $EXPERIMENT $sampleID"
    echo "$fullCmd"
    eval "$fullCmd"
    echo ""

    CURRENTSAMPLENUMBER=$((CURRENTSAMPLENUMBER+1))
done
echo "Done with experiment $EXPERIMENT"
CURRENTEXPNUMBER=$((CURRENTEXPNUMBER+1))
echo "Done with all experiments"