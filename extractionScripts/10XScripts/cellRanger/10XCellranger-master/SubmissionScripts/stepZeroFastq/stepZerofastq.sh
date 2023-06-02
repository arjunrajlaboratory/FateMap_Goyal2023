#!/bin/bash

# Change 
# input location is the location of the raw data - update this to be the location of your data in projects. fastqId is the folder which will be created after completion of this step and will contain the fastQ files.
# The final input is the csv file created which is specific to each experiment. Save this file in the repo folder of the project you are working on. 

projectName="projectABC"
experimentName="experimentXYZRun1"
fastqId="run1fastq"
csvName="SubmissionScripts/stepZeroFastq/exampleRun1.csv"

cellranger mkfastq --id=$fastqId --qc --localcores=16 \
                     --run=/project/arjunrajlab/$projectName/repo/10Xdatasets/$experimentName/bcl \
                     --csv=$csvName > /project/arjunrajlab/$projectName/repo/10Xdatasets/$experimentName/bcllog.txt