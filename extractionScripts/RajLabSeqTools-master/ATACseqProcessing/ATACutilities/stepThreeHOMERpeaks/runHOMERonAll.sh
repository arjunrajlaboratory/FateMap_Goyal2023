#!/bin/bash



# run from within repo 

# runs HOMER on each file

EXPERIMENT=$1

codeHomeDir=$2


CURRENTEXPNUMBER=1
CURRENTSAMPLENUMBER=1
for fileName in $EXPERIMENT/raw/* ; do
    echo "Working on sample $CURRENTSAMPLENUMBER of experiment $CURRENTEXPNUMBER"
    sampleID=`basename "$fileName"`

	# 
    fullCmd="$codeHomeDir/rajlabseqtools/ATACseqProcessing/ATACutilities/stepThreeHOMERpeaks/runHOMER.sh $EXPERIMENT $sampleID $STARFLAGS"
    echo "$fullCmd"
    eval "$fullCmd"
    echo ""

    CURRENTSAMPLENUMBER=$((CURRENTSAMPLENUMBER+1))
done
echo "Done with experiment $EXPERIMENT"
CURRENTEXPNUMBER=$((CURRENTEXPNUMBER+1))
echo "Done with all experiments"