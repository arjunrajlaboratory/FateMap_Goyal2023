#Script to parse timeMachine 100bp barcodes from single end FASTQ files off NextSeq.
#Requires the following python packages: Biopython, regex.

#Last updated 11/14/2018 to run on PMACS cluster

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip
import regex as re
from collections import Counter
from argparse import ArgumentParser
import os, glob
import numpy as np
from extractionFun import *

#Command line parser
parser = ArgumentParser()
#parser.add_argument("experiment", help = "Specify the path to the experiment directory")
parser.add_argument("sampleName", help = "Specify the name of sample directory containing the fastq.gz files")
parser.add_argument("-o", "--outFilePrefix", help = "Specify the output file prefix for table of barcodes and UMIs. If none specified, will use sampleDirectory") 
parser.add_argument("-s", "--stagger", help = "Specify the length of the stagger.", type=int, default = 0)
parser.add_argument("-r", "--includeReads", help = "If specified, output additional tables with barcodes and only read counts or UMI counts. Otherwise outputs only one table with both counts", action = 'store_true')
parser.add_argument("--check_vector", help = "Option to check vector sequence before or on both sides of the barcode sequence.", default = "both", choices = ["both", "before"]) 
parser.add_argument("-l", "--length", help = "If check_vector before specified, input here your desired barcode length. Default is 100", default = 100, type = int) 
parser.add_argument("-Q", "--minPhred", help = "Specify the minimum phredscore required to include a readout. Filters reads with more than 5 bases before the barcode with low phredscore.", default = 14, type = int) 
parser.add_argument("-e", "--excludedReads", help = "If specified, output txt.gz files containing reads excluded from the UMI and count files.", action = 'store_true')
args = parser.parse_args()

experimentDirectory = os.getcwd()

#Make directory for output files.
outFileDirectory = os.path.join(experimentDirectory, "analyzed", args.sampleName, 'extractedBarcodeData') 
if not os.path.exists(outFileDirectory):
	os.makedirs(outFileDirectory)

if args.outFilePrefix is not None:
	outFilePrefix = args.outFilePrefix
else:
	outFilePrefix = args.sampleName

if args.check_vector == 'both':
	outFileUMI = outFilePrefix + "_UMIs.gz"
	outFileCounts = outFilePrefix + "_counts.gz"
	outFileUMICounts = outFilePrefix + "_UMICountsOnly.gz"
	outFileReadCounts = outFilePrefix + "_readCountsOnly.gz"
elif args.check_vector == 'before':
	outFileUMI = outFilePrefix + "_UMIs_liberal.gz"
	outFileCounts = outFilePrefix + "_counts_liberal.gz"
	outFileUMICounts = outFilePrefix + "_UMICountsOnly_liberal.gz"
	outFileReadCounts = outFilePrefix + "_readCountsOnly_liberal.gz"

outFileMissingBeforeBarcode = outFilePrefix + "_missingBeforeBarcode.gz"
outFileMissingAfterBarcode = outFilePrefix + "_missingAfterBarcode.gz"
outFileBadLength = outFilePrefix + "_badLength.gz"
outFileBadPhred = outFilePrefix + "_badPhred.gz"
staggerLength = args.stagger
minPhred = args.minPhred

#Move to sample directory 
os.chdir(os.path.join(experimentDirectory, "raw", args.sampleName))

print "Parsing sample {}".format(args.sampleName)
inFileNames = glob.glob("*fastq*")
if args.check_vector == "both":
	barcode_dict, missingBeforeBarcode, missingAfterBarcode, badQscore, badLength = parseBarcodeFastQ_both(inFileNames, args.stagger, minPhred)
elif args.check_vector == "before":
	barcode_dict, missingBeforeBarcode, badQscore, badLength = parseBarcodeFastQ_before(inFileNames, args.stagger, args.length, minPhred)
print "Finished parsing sample {}. Length of dictionary is {}".format(args.sampleName, len(barcode_dict))  
print "Number of reads missing sequence before barcode is {}".format(len(missingBeforeBarcode))

if args.check_vector == "both":
	print "Number of reads missing sequence after barcode is {}".format(len(missingAfterBarcode))



#Write out barcode and associated phredscore and UMIs to file. 
os.chdir(outFileDirectory)
writeOutFileUMIs(barcode_dict, outFileUMI)             
print "Finished writing barcode UMIs to file"

if args.excludedReads == True:
	print "Writing exlcuded reads to files"
	writeOutFileBadSeqRecord(missingBeforeBarcode, outFileMissingBeforeBarcode)
	writeOutFileBadSeqRecord(badQscore, outFileBadPhred)
	writeOutFileBadSeqRecord(badLength, outFileBadLength)
	if args.check_vector == "both":
		writeOutFileBadSeqRecord(missingAfterBarcode, outFileMissingAfterBarcode)
    
print "Counting barcode UMIs"      
barcode_counts_dict = countUMIs(barcode_dict)                  
print "Finishing counting barcode UMIs. Length of dictionary:"  
print str(len(barcode_counts_dict))
    
#Write out barcode and associated read/UMI counts to file. 
if args.includeReads:
	writeOutFileBarcodeUMICounts(barcode_counts_dict,outFileUMICounts)
	writeOutFileBarcodeCounts(barcode_counts_dict, outFileCounts)
	writeOutFileBarcodeReadCounts(barcode_counts_dict, outFileReadCounts)
else:
	writeOutFileBarcodeCounts(barcode_counts_dict, outFileCounts)
    