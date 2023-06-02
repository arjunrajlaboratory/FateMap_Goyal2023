#!/bin/bash

# Change 
# input location is the location of the fast files.
projectName="projectABC"
experimentName="experimentXYZRun1"

sample1="/project/arjunrajlab/$projectName/repo/10Xdatasets/$experimentName/SubmissionScripts/stepOneCounts/stepOnecounts_s1.sh"
sample2="/project/arjunrajlab/$projectName/repo/10Xdatasets/$experimentName/SubmissionScripts/stepOneCounts/stepOnecounts_s2.sh"


JOURNAL="/project/arjunrajlab/$projectName/repo/10Xdatasets/$experimentName/submitStepOne_$(date +%Y-%m-%d_%H-%M).log"

bsub -n 16 -M 150000 -R "span[hosts=1] rusage [mem=150000]" < $sample1
bsub -n 16 -M 150000 -R "span[hosts=1] rusage [mem=150000]" < $sample2

