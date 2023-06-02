#!/bin/bash

# navigate to the directory right above the experiment directory
# run this script from there
EXPERIMENT="ProcessedDataEGFRsort"
#EXPERIMENT="ProcessedDataEGFRsort/differentialPeakCalls"

SAMPLE1="EGFR-week4"
SAMPLE2="EGFR-week1"

# R1 and R2 are replicates for each sample.
SAMPLE1_R1="EGFR-p9-week4"
SAMPLE1_R2="EGFR-p14-week4-v2"
SAMPLE2_R1="EGFR-p9-week1"
SAMPLE2_R2="EGFR-p14-week1-v2"

foldChange=4  # we ran this as foldChange 4 and 10
pval=0.0001 #this is default

# build out directory structure as needed...

#directory names:
DIFF="diffPeaks"
MERGE="mergePeaks"
ANNOTATE="annotatePeaks"
MOTIFS="motifs"

#PROCESSING="$SAMPLE1\Vs$SAMPLE2"
PROCESSING=$SAMPLE1\Vs$SAMPLE2

if [ ! -d $EXPERIMENT/$PROCESSING ]; then
    mkdir $EXPERIMENT/$PROCESSING/
    echo $EXPERIMENT/$PROCESSING/
fi


if [ ! -d $EXPERIMENT/$PROCESSING/foldChange$foldChange/ ]; then
    mkdir $EXPERIMENT/$PROCESSING/foldChange$foldChange/
    mkdir $EXPERIMENT/$PROCESSING/foldChange$foldChange/$DIFF
    mkdir $EXPERIMENT/$PROCESSING/foldChange$foldChange/$MERGE
    mkdir $EXPERIMENT/$PROCESSING/foldChange$foldChange/$ANNOTATE
    mkdir $EXPERIMENT/$PROCESSING/foldChange$foldChange/$MOTIFS
fi

# differential peak calling for R1:
getDifferentialPeaks $EXPERIMENT/analyzed/$SAMPLE1_R1/HOMER/peaks.txt \
$EXPERIMENT/analyzed/$SAMPLE1_R1/HOMER \
$EXPERIMENT/analyzed/$SAMPLE2_R1/HOMER \
-F $foldChange -P $pval -size 200 > \
$EXPERIMENT/$PROCESSING/foldChange$foldChange/$DIFF/diffPeaksR1_UpIn$SAMPLE1_R1.txt

# differential peak calling for R2:
getDifferentialPeaks $EXPERIMENT/analyzed/$SAMPLE1_R2/HOMER/peaks.txt \
$EXPERIMENT/analyzed/$SAMPLE1_R2/HOMER \
$EXPERIMENT/analyzed/$SAMPLE2_R2/HOMER \
-F $foldChange -P $pval -size 200 > \
$EXPERIMENT/$PROCESSING/foldChange$foldChange/$DIFF/diffPeaksR2_UpIn$SAMPLE1_R2.txt

## merge those 2 together...
cd $EXPERIMENT/$PROCESSING/foldChange$foldChange/
mergePeaks -prefix $MERGE/ \
-d given \
-venn $MERGE/Merge-UpIn$SAMPLE1-Venn.txt \
$DIFF/diffPeaksR1_UpIn$SAMPLE1_R1.txt \
$DIFF/diffPeaksR2_UpIn$SAMPLE1_R2.txt
cd ..
cd ..
cd ..

# _diffPeaks_diffPeaksR1-EGFR-week4Vsmix-noDrug.txt_diffPeaks_diffPeaksR2-EGFR-week4Vsmix-noDrug.txt

# getDifferentialPeaks only runs one direction, so to test the other direction, meaning
# find peaks that are present in sample 2, but not in sample 1, we must now run everything 
# swapped.

# differential peak calling for R1:
getDifferentialPeaks $EXPERIMENT/analyzed/$SAMPLE2_R1/HOMER/peaks.txt \
$EXPERIMENT/analyzed/$SAMPLE2_R1/HOMER \
$EXPERIMENT/analyzed/$SAMPLE1_R1/HOMER \
-F $foldChange -P $pval -size 200 > \
$EXPERIMENT/$PROCESSING/foldChange$foldChange/$DIFF/diffPeaksR1_UpIn$SAMPLE2_R1.txt

# differential peak calling for R2:
getDifferentialPeaks $EXPERIMENT/analyzed/$SAMPLE2_R2/HOMER/peaks.txt \
$EXPERIMENT/analyzed/$SAMPLE2_R2/HOMER \
$EXPERIMENT/analyzed/$SAMPLE1_R2/HOMER \
-F $foldChange -P $pval -size 200 > \
$EXPERIMENT/$PROCESSING/foldChange$foldChange/$DIFF/diffPeaksR2_UpIn$SAMPLE2_R2.txt

## merge those 2 together...
cd $EXPERIMENT/$PROCESSING/foldChange$foldChange/
mergePeaks -prefix $MERGE/ \
-d given \
-venn $MERGE/Merge-UpIn$SAMPLE2-Venn.txt \
$DIFF/diffPeaksR1_UpIn$SAMPLE2_R1.txt \
$DIFF/diffPeaksR2_UpIn$SAMPLE2_R2.txt
cd ..
cd ..
cd ..


# rename Files....
# _diffPeaks_diffPeaksR1_UpIn$SAMPLE1_R1.txt_diffPeaks_diffPeaksR2_UpInE$SAMPLE1_R2.txt
mv $EXPERIMENT/$PROCESSING/foldChange$foldChange/$MERGE/_diffPeaks_diffPeaksR1_UpIn$SAMPLE1_R1.txt_diffPeaks_diffPeaksR2_UpIn$SAMPLE1_R2.txt \
$EXPERIMENT/$PROCESSING/foldChange$foldChange/$MERGE/Merge-UpIn$SAMPLE1.txt

mv $EXPERIMENT/$PROCESSING/foldChange$foldChange/$MERGE/_diffPeaks_diffPeaksR1_UpIn$SAMPLE2_R1.txt_diffPeaks_diffPeaksR2_UpIn$SAMPLE2_R2.txt \
$EXPERIMENT/$PROCESSING/foldChange$foldChange/$MERGE/Merge-UpIn$SAMPLE2.txt


# annotate files...
annotatePeaks.pl $EXPERIMENT/$PROCESSING/foldChange$foldChange/$MERGE/Merge-UpIn$SAMPLE1.txt \
hg19 -go $EXPERIMENT/$PROCESSING/foldChange$foldChange/$ANNOTATE/GeneOntology_UpIn$SAMPLE1/ \
-genomeOntology $EXPERIMENT/$PROCESSING/foldChange$foldChange/$ANNOTATE/GenomeOntology_UpIn$SAMPLE1/ \
> $EXPERIMENT/$PROCESSING/foldChange$foldChange/$ANNOTATE/Annotate-UpIn$SAMPLE1.txt

annotatePeaks.pl $EXPERIMENT/$PROCESSING/foldChange$foldChange/$MERGE/Merge-UpIn$SAMPLE2.txt \
hg19 -go $EXPERIMENT/$PROCESSING/foldChange$foldChange/$ANNOTATE/GeneOntology_UpIn$SAMPLE2/ \
-genomeOntology $EXPERIMENT/$PROCESSING/foldChange$foldChange/$ANNOTATE/GenomeOntology_UpIn$SAMPLE2/ \
> $EXPERIMENT/$PROCESSING/foldChange$foldChange/$ANNOTATE/Annotate-UpIn$SAMPLE2.txt


# find motifs - up in sample 1
findMotifsGenome.pl \
$EXPERIMENT/$PROCESSING/foldChange$foldChange/$ANNOTATE/Annotate-UpIn$SAMPLE1.txt \
hg19 $EXPERIMENT/$PROCESSING/foldChange$foldChange/$MOTIFS/Annotate-UpIn$SAMPLE1.txt \
-size given

# find motifs - up in sample 2
findMotifsGenome.pl \
$EXPERIMENT/$PROCESSING/foldChange$foldChange/$ANNOTATE/Annotate-UpIn$SAMPLE2.txt \
hg19 $EXPERIMENT/$PROCESSING/foldChange$foldChange/$MOTIFS/Annotate-UpIn$SAMPLE2.txt \
-size given


	