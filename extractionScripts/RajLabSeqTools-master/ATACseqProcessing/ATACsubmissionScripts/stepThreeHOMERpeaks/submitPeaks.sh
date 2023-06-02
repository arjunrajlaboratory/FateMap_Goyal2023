#!/bin/bash


EXPERIMENT="ProcessedDataEGFRsort"

codeHomeDir="/Users/sydneyshaffer"

cmdToRun="$codeHomeDir/rajlabseqtools/ATACseqProcessing/ATACutilities/stepThreeHOMERpeaks/runHOMERonAll.sh \
	$EXPERIMENT $codeHomeDir"

echo "$cmdToRun"
eval "$cmdToRun"
echo "done"