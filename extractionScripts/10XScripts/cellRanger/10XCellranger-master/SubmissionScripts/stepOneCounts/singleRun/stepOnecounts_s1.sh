#!/bin/bash

# Change 
# input location is the location of the fast files.
# sampleID is same as what was input in the CSV file in stepZero

projectName="projectABC"
experimentName="experimentXYZRun1"
sampleID=run1_sample1
fastqR="/project/arjunrajlab/$projectName/repo/10Xdatasets/$experimentName/run1fastq/outs/fastq_path"
refGenome="/home/goyaly/code/refdata-gex-GRCh38-2020-A"
log="countlogr1s1.txt"

JOURNAL="/project/arjunrajlab/$projectName/repo/10Xdatasets/$experimentName/$(date +%Y-%m-%d_%H-%M).log"

	
	cellranger count --id=$sampleID-count --localcores=16 \
                   --transcriptome=$refGenome \
                   --fastqs=$fastqR \
                   --sample=$sampleID \
                   --expect-cells=9000 > /project/arjunrajlab/$projectName/repo/10Xdatasets/$experimentName/outlogs/$log


	    	date >> $JOURNAL
    	echo "Done" >> $JOURNAL
