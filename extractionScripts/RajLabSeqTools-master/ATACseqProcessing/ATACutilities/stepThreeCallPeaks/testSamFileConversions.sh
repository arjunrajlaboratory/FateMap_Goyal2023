#!/bin/bash


EXPERIMENT="ProcessedDataEGFRsort"

sampleID="EGFR-p14-noDrug-v2"

codeHomeDir="/Users/sydneyshaffer"

fullCmd="$codeHomeDir/rajlabseqtools/ATACseqProcessing/ATACutilities/stepThreeCallPeaks/samToSortedBam.sh $EXPERIMENT $sampleID"

echo "$fullCmd"
eval "$fullCmd"