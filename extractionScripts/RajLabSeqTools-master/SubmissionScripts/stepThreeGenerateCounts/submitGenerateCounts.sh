#!
# run from within repo

cmdToRun="$codeHomeDir/rajlabseqtools/Utilities/stepThreeGenerateCounts/allGenerateCounts.sh $EXPERIMENT $codeHomeDir $PROJECT $gtfFile"

JOB_ARRAY_PARAMETER_TEMPLATE="generateCounts[1-tk]"

#The expression within ${} below replaces "tk" with the correct number of samples (defined in the setEnvironmentVariables.sh script)
bsub -J ${JOB_ARRAY_PARAMETER_TEMPLATE/tk/$N_SAMPLES} -o out.%I -e err.%I $cmdToRun
