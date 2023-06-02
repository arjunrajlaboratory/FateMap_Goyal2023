#!/bin/bash
# run from within repo 

cmdToRun="$codeHomeDir/rajlabseqtools/Utilities/stepTwoStar/runStarOnAll.sh $EXPERIMENT $codeHomeDir $PAIRED_OR_SINGLE_END_FRAGMENTS $genomeDirSTAR $STARFLAGS"

echo "$cmdToRun"
eval "$cmdToRun"
echo "done submitting star jobs, wait until all jobs are complete before proceeding to the next step. use the bjobs command to monitor their progress."

