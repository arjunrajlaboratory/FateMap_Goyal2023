#!/bin/bash
#BSUB -J bowtie2
#BSUB -o bowtie2.%J.out
#BSUB -e bowtie2.%J.error



# run from within repo 

EXPERIMENT=$1

codeHomeDir=$2

STARFLAGS=${@:3} # pass all arguments after the first two

CURRENTEXPNUMBER=1
CURRENTSAMPLENUMBER=1
for fileName in $EXPERIMENT/raw/* ; do
    echo "Working on sample $CURRENTSAMPLENUMBER of experiment $CURRENTEXPNUMBER"
    sampleID=`basename "$fileName"`

    fullCmd="$codeHomeDir/rajlabseqtools/ATACseqProcessing/ATACutilities/stepTwoBowtie2/runBowtie2.sh $EXPERIMENT $sampleID $STARFLAGS"
    echo "$fullCmd"
    eval "$fullCmd"
    echo ""

    CURRENTSAMPLENUMBER=$((CURRENTSAMPLENUMBER+1))
done
echo "Done with experiment $EXPERIMENT"
CURRENTEXPNUMBER=$((CURRENTEXPNUMBER+1))
echo "Done with all experiments"



#echo "Renaming SAM files to match sample names"
#fullCmd2="$codeHomeDir/rajlabseqtools/ATACseqProcessing/Utilities/stepTwoStar/makeSamRenamingLinks.sh $EXPERIMENT"
#echo "$fullCmd2"
#eval "$fullCmd2"
#echo "done renaming"
