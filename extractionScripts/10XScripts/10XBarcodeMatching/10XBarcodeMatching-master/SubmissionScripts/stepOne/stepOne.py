#==================================================================================================
#Input Files: This pipeline takes fastq.gz files that Illumina returns as input. 
#Change these file PATHs based on your folder structure and where your datasets are stored. 
#$PATH needs to be changed at 3 places; Stagger depends on the i7 primer and must be specified.
#==================================================================================================

inputFastqGzReadR1 = "/project/arjunrajlab/experimentName/Analysis/stepZero/raw/s1/s1.fastq.gz"
inputFastqGzReadR2 = "/project/arjunrajlab/experimentName/Analysis/stepZero/raw/s1/s1.fastq.gz"
staggerLength = 1 #length of the stagger (e.g. For i7.3.1 corresponding to sample 1, it is 1)
outputDirectory = "/project/arjunrajlab/experimentName/Analysis/stepOne/s1/"

#=================================================================================================
#*************************************************************************************************
#******************************DO NOT EDIT ANYTHING BELOW THIS POINT *****************************
#*************************************************************************************************		
#=================================================================================================

from Bio import SeqIO
import gzip
import os, sys
import subprocess
import regex
from gzip import open as gzopen
import numpy as np

read1 = []  #read 1, corresponding to CellID
read2 = []  #read 2, corresponding to GFP+ barcode
read1Qscore = [] #read 1 quality scores as array
read2Qscore = [] #read 2 quality scores as array
joinedRead1Read2 =[]
missingGfpBarcode = []
badQscore = []
screenedReads = []
shavedReads = []
badBarcode = []
testConditionID = []
uniqueScreenedReads = []
uniqueShavedReads = []
gfpPrimerLength = 20 #including (update - excluding) two extra bases AT after the GFP primer
recordsR1 = SeqIO.parse(gzopen(inputFastqGzReadR1, "rt"),format="fastq")
recordsR2 = SeqIO.parse(gzopen(inputFastqGzReadR2, "rt"),format="fastq")
sumGfpPrimerStaggrer = gfpPrimerLength + staggerLength

print("Matching reads...")
# Automatic check on whether the identifiers match for Read 1 and Read 2     
for record1, record2 in zip(recordsR1, recordsR2): 
        if record1.id == record2.id:
            read1.append(str(record1.seq))
            read2.append(str(record2.seq))
            read1Qscore.append(record1.letter_annotations["phred_quality"])
            read2Qscore.append(record2.letter_annotations["phred_quality"])

read2QscoreInt = np.array(read2Qscore) # convert it to array otherwise it gives error
print("finished matching reads")            
joinedRead1Read2 = np.column_stack(((read1,read2))) # joining read 1 and read 2
print("joinRead1Read2Done")

print("Starting filtering steps...")
for read, Qscore in zip(joinedRead1Read2, read2QscoreInt):
	read1Call = read[0]
	read2Call = read[1]
	cellID = read[0][0:16]
	UMI = read[0][17:26]
	barcode = read[1][sumGfpPrimerStaggrer:(sumGfpPrimerStaggrer + 70)]  #needs to change depending on the sample
	shaved = np.array([cellID, UMI,barcode])
	read2QscoreCall = Qscore
	conditionGFP = 1
	conditionQscore = 1
	conditionBadBarcode = 1

	if len(regex.findall(r'(GGACGAGCTGTACAAGTAGG){e<=4}', read2Call[0:(sumGfpPrimerStaggrer)])) == 0: # Do not consider reads without GFP tag
		missingGfpBarcode.append(read)
		conditionGFP = 0
	if len(read2QscoreCall[np.where(read2QscoreCall[0:(sumGfpPrimerStaggrer)] <= 14)]) > 5: # Do not consider reads with phred Q score <14 for more than 5 times in GFP
		badQscore.append(read)
		conditionQscore = 0
	if (len(regex.findall("(AAAA)", read2Call)) > 0 or \
        len(regex.findall("(TTTT)", read2Call)) > 0 or \
        len(regex.findall("(GGGG)", read2Call)) > 0 or \
        len(regex.findall("(CCCC)", read2Call)) > 0 or \
        len(regex.findall("(NN)", read2Call)) > 0):
		badBarcode.append(read)
		conditionBadBarcode = 0
	if (conditionQscore == 1 and conditionGFP == 1 and conditionBadBarcode == 1):
		screenedReads.append(read)
		shavedReads.append(shaved)
        
uniqueScreenedReads = np.unique(screenedReads, axis = 0)
uniqueShavedReads = np.unique(shavedReads, axis = 0)

summaryFile = open(outputDirectory + "summaryFile.txt","w") 
summaryFile.write("total raw reads %d" % len(read1) + "\n")
summaryFile.write("total missingGfpBarcode reads %d" % len(missingGfpBarcode) + "\n") 
summaryFile.write("total badQscore reads %d" % len(badQscore) + "\n") 
summaryFile.write("total badBarcode reads %d" % len(badBarcode) + "\n") 
summaryFile.write("total screenedReads reads %d" % len(screenedReads) + "\n")
summaryFile.write("total screened Unique Reads reads %d" % len(uniqueScreenedReads) + "\n") 
summaryFile.write("total shavedReads reads %d" % len(shavedReads) + "\n")
summaryFile.write("total shaved Unique Reads reads %d" % len(uniqueShavedReads) + "\n") 
summaryFile.close()


#writing files to the Analysis folder
np.savetxt(outputDirectory + "joinedRead1Read2.txt", joinedRead1Read2, fmt='%s')
np.savetxt(outputDirectory + "badQscore.txt", badQscore, fmt='%s')  
np.savetxt(outputDirectory + "badBarcode.txt", badBarcode, fmt='%s')  
np.savetxt(outputDirectory + "uniqueScreenedReads.txt", uniqueScreenedReads, fmt='%s')  
np.savetxt(outputDirectory + "uniqueShavedReads.txt", uniqueShavedReads, fmt='%s')  



