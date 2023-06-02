#!/bin/bash

ZIPFILEDIRECTORY=$1

OUTFASTQDIRECTORY=$2

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

    FASTQR1=${SAMPLE}.fastq
#    FASTQR2=${SAMPLE}_R2.fastq

    if [ ! -e $OUTFASTQDIRECTORY/raw/$SAMPLE/$FASTQR1 ]; then
        echo Working on $SAMPLE
        for i in *.gz; do
            gunzip -c $i > ${i%.*}
        done


        cat *fastq > $OUTFASTQDIRECTORY/raw/$SAMPLE/$FASTQR1
    #    cat *R2*fastq > $OUTFASTQDIRECTORY/raw/$SAMPLE/$FASTQR2
    fi
done

