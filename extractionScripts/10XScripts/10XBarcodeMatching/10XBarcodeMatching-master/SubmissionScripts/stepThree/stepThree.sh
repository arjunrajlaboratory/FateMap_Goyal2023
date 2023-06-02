#!/bin/bash

# input location is the location of the raw data - update this to be the location of your data in projects. fastqId is the folder which will be created after completion of this step and will contain the fastQ files.
# The final input is the csv file created which is specific to each experiment. Save this file in the repo folder of the project you are working on. 
# LV distance can be specified (-d). Default is 8 (-d8) for three files, 6 for one file. 
inputDirectory="/project/arjunrajlab/experimentName/Analysis/stepTwo"
outputDirectory="/project/arjunrajlab/experimentName/Analysis/stepThree"
PATH=$PATH:/home/goyaly/code/starcode

printf "starcode running\n"

starcode -i $inputDirectory/stepTwoBarcodes50.txt -d8 -o $outputDirectory/stepThreeBarcodes50_d8 --seq-id -s > $outputDirectory/50_8log.txt
printf "50_d8 done\n"
starcode -i $inputDirectory/stepTwoBarcodes40.txt -d8 -o $outputDirectory/stepThreeBarcodes40_d8 --seq-id -s > $outputDirectory/40_8log.txt
printf "40_d8 done\n"
starcode -i $inputDirectory/stepTwoBarcodes30.txt -d8 -o $outputDirectory/stepThreeBarcodes30_d8 --seq-id -s > $outputDirectory/30_8log.txt
printf "30_d8 done\n"
starcode -i $inputDirectory/stepTwoBarcodes30.txt -d6 -o $outputDirectory/stepThreeBarcodes30_d6 --seq-id -s > $outputDirectory/30_6log.txt
printf "30_d6 done\n"

python stepThreeRun.py $inputDirectory/ $outputDirectory/