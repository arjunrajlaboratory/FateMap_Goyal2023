#!/bin/bash

ZIPFILEDIRECTORY=$1
OUTFASTQDIRECTORY=$2
PAIRED_OR_SINGLE_END_FRAGMENTS=$3

echo "you selected the pipeline settings for" $PAIRED_OR_SINGLE_END_FRAGMENTS "end reads" 

for dirname in $ZIPFILEDIRECTORY/* ; do
    cd $dirname

    INPUT=`ls *`
    SAMPLE=${INPUT%%_*}  # Cuts filename string after first '_'

    if [ ! -d $OUTFASTQDIRECTORY/raw ]; then
        mkdir $OUTFASTQDIRECTORY/raw
    fi

    if [ ! -d $OUTFASTQDIRECTORY/raw/$SAMPLE ]; then
        mkdir $OUTFASTQDIRECTORY/raw/$SAMPLE
    fi

    FASTQR1=${SAMPLE}_R1.fastq
    FASTQR2=${SAMPLE}_R2.fastq

    if [ ! -e $OUTFASTQDIRECTORY/raw/$SAMPLE/$FASTQR1 ]; then
        echo Working on $SAMPLE
        for i in *.gz; do
            gunzip -c $i > ${i%.*}
        done


        cat ./*R1*fastq > $OUTFASTQDIRECTORY/raw/$SAMPLE/$FASTQR1
        rm ./*R1*fastq 

        if [ $PAIRED_OR_SINGLE_END_FRAGMENTS = "paired" ]; then
            cat ./*R2*fastq > $OUTFASTQDIRECTORY/raw/$SAMPLE/$FASTQR2
            rm ./*R2*fastq
        fi
    fi
    
    # compress fastq after concatenation, since STAR can read compressed FASTQ files.
    # submit this final compression task to cluster since compression takes a while
    bsub gzip $OUTFASTQDIRECTORY/raw/$SAMPLE/$FASTQR1
    if [ $PAIRED_OR_SINGLE_END_FRAGMENTS = "paired" ]; then
        bsub gzip $OUTFASTQDIRECTORY/raw/$SAMPLE/$FASTQR2
    fi
done

