#!/bin/bash

# Change 
# input location is the location of the fast files.

projectName="projectABC"
experimentName="experimentXYZRun1Run2"
experimentName1="experimentXYZRun1"
experimentName2="experimentXYZRun2"
sampleID1=run1_sample2
sampleID2=run2_sample2
fastqR1="/project/arjunrajlab/$projectName/repo/10Xdatasets/$experimentName1/run1fastq/outs/fastq_path"
fastqR2="/project/arjunrajlab/$projectName/repo/10Xdatasets/$experimentName2/run2fastq/outs/fastq_path"
refGenome="/home/goyaly/code/refdata-gex-GRCh38-2020-A"
log="countlogr12s2.txt"

JOURNAL="/project/arjunrajlab/$projectName/repo/10Xdatasets/$experimentName/$(date +%Y-%m-%d_%H-%M).log"

	
	cellranger count --id=$sampleID1-$sampleID2-count --localcores=16 \
                   --transcriptome=$refGenome \
                   --fastqs=$fastqR1,$fastqR2 \
                   --sample=$sampleID1,$sampleID2 \
                   --expect-cells=9000 > /project/arjunrajlab/$projectName/repo/10Xdatasets/$experimentName/outlogs/$log


	    	date >> $JOURNAL
    	echo "Done" >> $JOURNAL
