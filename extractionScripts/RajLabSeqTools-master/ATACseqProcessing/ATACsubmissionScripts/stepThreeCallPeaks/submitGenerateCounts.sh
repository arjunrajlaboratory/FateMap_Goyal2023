#!
# run from within repo
projectName="CancerSeq"

scriptFile=/project/arjunrajlab/$projectName/repo/SubmissionScripts/stepThreeGenerateCounts/defineVariablesGenerateCounts.sh
bsub -J "generateCounts[1-24]" -o out.%I -e err.%I < $scriptFile
