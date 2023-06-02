#!/bin/bash


EXPERIMENT="ProcessedDataEGFRsort"

codeHomeDir="/Users/sydneyshaffer"

cmdToRun="$codeHomeDir/rajlabseqtools/ATACseqProcessing/ATACutilities/stepThreeCallPeaks/runPeakCallersOnAll.sh \
	$EXPERIMENT $codeHomeDir"

echo "$cmdToRun"
eval "$cmdToRun"
echo "done"