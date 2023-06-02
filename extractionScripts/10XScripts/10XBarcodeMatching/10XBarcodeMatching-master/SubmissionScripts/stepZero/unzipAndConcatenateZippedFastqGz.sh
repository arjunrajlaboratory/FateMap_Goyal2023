#!/bin/bash

ZIPFILEDIRECTORY=$1

OUTFASTQDIRECTORY=$2

for dirname in $ZIPFILEDIRECTORY/* ; do
    cd $dirname

    INPUT=`ls *`
    SAMPLE=${INPUT%%_*}  # Cuts filename string after first '_'

    if [ ! -d $OUTFASTQDIRECTORY/$SAMPLE ]; then
        mkdir $OUTFASTQDIRECTORY/$SAMPLE
    fi

    FASTQR1=${SAMPLE}_R1.fastq
    FASTQR2=${SAMPLE}_R2.fastq

    if [ ! -e $OUTFASTQDIRECTORY/$SAMPLE/$FASTQR1 ]; then
        echo Working on $SAMPLE
        for i in *.gz; do
            gunzip -c $i > ${i%.*}
        done

        cat *R1*fastq > $OUTFASTQDIRECTORY/$SAMPLE/$FASTQR1
        cat *R2*fastq > $OUTFASTQDIRECTORY/$SAMPLE/$FASTQR2

	gzip $OUTFASTQDIRECTORY/$SAMPLE/$FASTQR1 	#to zip again for this pipeline
	gzip $OUTFASTQDIRECTORY/$SAMPLE/$FASTQR2	#to zip again for this pipeline
    fi
done
