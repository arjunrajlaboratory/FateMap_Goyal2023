#!/bin/bash


EXPERIMENT="ProcessedDataEGFRsort"
STARFLAGS=""

codeHomeDir=/Users/sydneyshaffer

cmdToRun="$codeHomeDir/rajlabseqtools/ATACseqProcessing/ATACutilities/stepTwoBowtie2/runBowtie2OnAll.sh $EXPERIMENT $codeHomeDir $STARFLAGS"

echo "$cmdToRun"
eval "$cmdToRun"
echo "done with bowtie2"

