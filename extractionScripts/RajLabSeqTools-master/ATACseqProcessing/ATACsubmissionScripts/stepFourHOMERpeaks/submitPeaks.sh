#!/bin/bash


EXPERIMENT="ProcessedDataEGFRsort"

codeHomeDir="~"

cmdToRun="$codeHomeDir/rajlabseqtools/ATACseqProcessing/ATACutilities/stepFourHOMERpeaks/runHOMERonAll.sh \
	$EXPERIMENT $codeHomeDir"

echo "$cmdToRun"
eval "$cmdToRun"
echo "done"