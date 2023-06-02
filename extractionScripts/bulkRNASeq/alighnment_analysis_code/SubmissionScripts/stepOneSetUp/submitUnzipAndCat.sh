#!/bin/bash

# The first input is the location of the raw data - this should be the location of your data.
# The second/final input is the destination - this should be your experiment folder.

$codeHomeDir/rajlabseqtools/Utilities/stepOneSetUp/unzipAndConcatenateZippedFastq.sh $RAWDATA_DIRECTORY /project/arjunrajlab/$PROJECT/repo/$EXPERIMENT $PAIRED_OR_SINGLE_END_FRAGMENTS