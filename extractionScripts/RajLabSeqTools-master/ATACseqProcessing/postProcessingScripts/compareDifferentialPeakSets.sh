#!/bin/bash

# we want to know how many of the differential peaks between week 1 and week 4 are the 
# same differential peaks that we saw between non-resistant and week 4

EXPERIMENT="ProcessedDataEGFRsort/differentialPeakCalls"

PEAKS_W1_W4="$EXPERIMENT/EGFR-week4VsEGFR-week1/foldChange4/mergePeaks/Merge-UpInEGFR-week4.txt"

PEAKS_nonR_W4="$EXPERIMENT/EGFR-week4Vsmix-noDrug/foldChange4/mergePeaks/Merge-UpInEGFR-week4.txt"

pathToOutput="$EXPERIMENT/postProcessing/merge_Week1VsWeek4_Week4VsNonR"

mergePeaks -prefix $pathToOutput/ \
-d 1000 \
-venn $pathToOutput/Venn.txt \
$PEAKS_W1_W4 $PEAKS_nonR_W4 

# rename these files
mv $pathToOutput/_ProcessedDataEGFRsort_differentialPeakCalls_EGFR-week4VsEGFR-week1_foldChange4_mergePeaks_Merge-UpInEGFR-week4.txt_ProcessedDataEGFRsort_differentialPeakCalls_EGFR-week4Vsmix-noDrug_foldChange4_mergePeaks_Merge-UpInEGFR-week4.txt \
$pathToOutput/differentialInBoth.txt

mv $pathToOutput/_ProcessedDataEGFRsort_differentialPeakCalls_EGFR-week4Vsmix-noDrug_foldChange4_mergePeaks_Merge-UpInEGFR-week4.txt \
$pathToOutput/differentialInWeek4VsNonR.txt

mv $pathToOutput/_ProcessedDataEGFRsort_differentialPeakCalls_EGFR-week4VsEGFR-week1_foldChange4_mergePeaks_Merge-UpInEGFR-week4.txt \
$pathToOutput/differentialInWeek4VsWeek1.txt


#############

# now we want to make the same type of comparison, but for the peaks that go away
# here it is peaks nonR vs week 1 and nonR vs week4

PEAKS_nonR_W1="$EXPERIMENT/EGFR-week1Vsmix-noDrug/foldChange4/mergePeaks/Merge-UpInmix-noDrug.txt"

PEAKS_nonR_W4="$EXPERIMENT/EGFR-week4Vsmix-noDrug/foldChange4/mergePeaks/Merge-UpInmix-noDrug.txt"

pathToOutput="$EXPERIMENT/postProcessing/merge_Week1VsNonR_Week4VsNonR"

mergePeaks -prefix $pathToOutput/ \
-d 1000 \
-venn $pathToOutput/Venn.txt \
$PEAKS_nonR_W1 $PEAKS_nonR_W4 

# rename these files
mv $pathToOutput/_ProcessedDataEGFRsort_differentialPeakCalls_EGFR-week1Vsmix-noDrug_foldChange4_mergePeaks_Merge-UpInmix-noDrug.txt_ProcessedDataEGFRsort_differentialPeakCalls_EGFR-week4Vsmix-noDrug_foldChange4_mergePeaks_Merge-UpInmix-noDrug.txt \
$pathToOutput/differentialInBoth.txt

mv $pathToOutput/_ProcessedDataEGFRsort_differentialPeakCalls_EGFR-week1Vsmix-noDrug_foldChange4_mergePeaks_Merge-UpInmix-noDrug.txt \
$pathToOutput/differentialInWeek1VsNonR.txt

mv $pathToOutput/_ProcessedDataEGFRsort_differentialPeakCalls_EGFR-week4Vsmix-noDrug_foldChange4_mergePeaks_Merge-UpInmix-noDrug.txt \
$pathToOutput/differentialInWeek4VsNonR.txt

###################

# next question: 
# are the peaks that go away at week 1 (but are not differential in nonR vs week4) the 
# same peaks that are coming back at week 4?

PEAKS_awayInWeek1="$EXPERIMENT/postProcessing/merge_Week1VsNonR_Week4VsNonR/differentialInWeek1VsNonR.txt"

PEAKS_upInWeek4="$EXPERIMENT/postProcessing/merge_Week1VsWeek4_Week4VsNonR/differentialInWeek4VsWeek1.txt"

pathToOutput="$EXPERIMENT/postProcessing/merge_PeaksAwayThenBack"

mergePeaks -prefix $pathToOutput/ \
-d 1000 \
-venn $pathToOutput/Venn.txt \
$PEAKS_awayInWeek1 $PEAKS_upInWeek4
