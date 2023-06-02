#!/bin/bash
# run in interactive session from within repo

echo "...generating melted data"
meltCmd="$codeHomeDir/rajlabseqtools/Utilities/stepThreeGenerateCounts/meltHTSeqData.pl $EXPERIMENT/analyzed/*/htseq/*.htseq.stdout > /project/arjunrajlab/$PROJECT/repo/$EXPERIMENT/${EXPERIMENT}_meltedData.tsv"
eval "$meltCmd"
echo "done"