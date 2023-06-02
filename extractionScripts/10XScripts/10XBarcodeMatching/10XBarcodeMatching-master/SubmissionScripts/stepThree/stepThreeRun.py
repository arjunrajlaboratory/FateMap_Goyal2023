
import os, sys
import subprocess
from gzip import open as gzopen
import numpy as np
import csv
import pandas as pd
    
csv.field_size_limit(100000000)

inputDirectory = sys.argv[1]
outputDirectory = sys.argv[2]
stepThreeStarcodeShavedReads = np.genfromtxt(inputDirectory + "stepTwoCellIDUMIBarcodes.txt", skip_header = 1, dtype = 'str')


print("Starcode files saved...")

with open(outputDirectory + "stepThreeBarcodes50_d8", 'r') as csvfile:
    stepThreeBarcodes50_d8 = csv.reader(csvfile, delimiter='\t')
    for i in stepThreeBarcodes50_d8:
        indices = [int(x) for x in i[2].split(',')]
        sequence = i[0]
        indices = np.array(indices)
        indices = indices-1
        stepThreeStarcodeShavedReads[indices,3] = sequence

with open(outputDirectory + "stepThreeBarcodes40_d8", 'r') as csvfile:
    stepThreeBarcodes40_d8 = csv.reader(csvfile, delimiter='\t')
    for i in stepThreeBarcodes40_d8:
        indices = [int(x) for x in i[2].split(',')]
        sequence = i[0]
        indices = np.array(indices)
        indices = indices-1
        stepThreeStarcodeShavedReads[indices,4] = sequence

with open(outputDirectory + "stepThreeBarcodes30_d8", 'r') as csvfile:
    stepThreeBarcodes30_d8 = csv.reader(csvfile, delimiter='\t')
    for i in stepThreeBarcodes30_d8:
        indices = [int(x) for x in i[2].split(',')]
        sequence = i[0]
        indices = np.array(indices)
        indices = indices-1
        stepThreeStarcodeShavedReads[indices,5] = sequence        
        
with open(outputDirectory + "stepThreeBarcodes30_d6", 'r') as csvfile:
    stepThreeBarcodes30_d6 = csv.reader(csvfile, delimiter='\t')
    for i in stepThreeBarcodes30_d6:
        indices = [int(x) for x in i[2].split(',')]
        sequence = i[0]
        indices = np.array(indices)
        indices = indices-1
        stepThreeStarcodeShavedReads[indices,6] = sequence               

print("Creating starcode replaced table...")

#writing files to the analysis folder, stepThree
colNames = ['cellID', 'UMI', 'OriginalBarcode', 'BC50StarcodeD8', 'BC40StarcodeD8', 'BC30StarcodeD8', 'BC30StarcodeD6', 'SampleNum']
stepThreeStarcodeShavedReads = pd.DataFrame(stepThreeStarcodeShavedReads, columns=colNames)
stepThreeStarcodeShavedReads.to_csv(outputDirectory + "stepThreeStarcodeShavedReads.txt", sep='\t', index=False)