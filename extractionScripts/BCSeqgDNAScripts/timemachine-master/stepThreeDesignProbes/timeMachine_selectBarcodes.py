#Python script for selecting time machine barcodes and writing barcodes to fasta files. 

import os, glob 
import itertools
import numpy as np
import pandas as pd
import Levenshtein
import regex as re
from argparse import ArgumentParser

#Command line parser
parser = ArgumentParser()
parser.add_argument("path", help = "Specify the path to the experiment directory containing subdirectories with the extracted barcodes.", type = str)
parser.add_argument("sample", help = "Specify the name(s) of the sample(s). Option to specify multiple replicate samples", nargs = '+', type = str)
parser.add_argument("number", help = "Specify the number of barcodes to write to fasta files for each sample", type = int)
parser.add_argument("-c", "--countMetric", help = "Specify weather to use read counts or 'UMI' counts to rank barcode sequence.", default = 'reads', choices = ['reads', 'UMIs'])
parser.add_argument("-j", "--join", help = "Option to join the samples (e.g. technical replicates", action = 'store_true')
parser.add_argument("-m", "--mask", help = "Option to specify a mask sequence (fasta format). Use absolute path. Will exclude barcodes that match up to numMask against the mask sequence. Input single fasta formatted sequence", type = str)
parser.add_argument("-n", "--numMask", help = "If mask sequence specified, specify the minimum hammning distance needed to mask", type = int)
parser.add_argument("-l", "--minLength", help = "Option to specify a minimum length for selected barcodes. Default is 90 bps", default = 90, type = int)
parser.add_argument("-o", "--outFileName", help = "Option to specify the prefix of the output files. Default will be to use the first input sample name", type = str)
parser.add_argument("--outputMasked", help = "Option to output barcodes that were masked to separate file", action = 'store_true')
parser.add_argument("-e", "--extend", help = "Option to add additional nucleotides to 5' and 3' ends of selected barcodes. Input 2 strings composed of A, C, G or T", nargs = 2, type = str)
parser.add_argument("-s", "--starcodeDistance", help = "Option to specify using barcode counts after running starcode with specified distance", type = int)
args = parser.parse_args()

def reverse_complement(sequence):
	complement = {'A':'T', 'T': 'A', 'G':'C', 'C':'G', 'N':'N'}
	reverseComplement = "".join([complement[base] for base in str(sequence)[::-1].upper()])
	return reverseComplement     

def writeToFasta(path, name, barcodeSeq):
	"""Function for pandas dataframe to write a barcode sequence to fasta file using barcode rank for 
	file name. BarcodeSeq should be a string of the entire barcode sequence. 
	"""
	with open(os.path.join(path, "barcode{}.fa".format(name+1)), 'w') as fasta: #Make barcode file and start at 1 rather than 0. 
		fasta.write(">barcode{}\n".format(name+1))
		fasta.write(barcodeSeq + "\n")

def getMinHamming(barcodeSequence, maskSequence):
	length = len(barcodeSequence)
	maskSequenceList = [maskSequence[i:i+length] for i in xrange(0, len(maskSequence)-length + 1)] 
	minHamming = min([Levenshtein.hamming(reverse_complement(barcodeSequence), j.upper()) for j in maskSequenceList])
	return minHamming

os.chdir(args.path)

samples = args.sample

if args.starcodeDistance is not None:
	countFiles = ['analyzed/{}/starcode/{}_clusteredBarcode{}_d{}.txt'.format(i, i, args.countMetric.capitalize(), args.starcodeDistance) for i in samples]

else:
	countFiles = ['analyzed/{}/extractedBarcodeData/{}_counts.gz'.format(i, i) for i in samples]

if args.starcodeDistance is not None:
	barcodeDataFrames = [pd.read_csv(i, names=["barcode", "counts"], sep = '\t') for i in countFiles]
else:
	barcodeDataFrames = [pd.read_csv(i, names=["barcode", "reads", "UMIs"], sep = '\t') for i in countFiles]

#Merge sample data frames if replicates. For simplicity, assumes there are only 2 samples to merge. Needs to be update to allow for more than 2 samples. 
if args.join == True:
	barcodeDf = pd.merge(barcodeDataFrames[0], barcodeDataFrames[1], on = 'barcode')
	if args.starcodeDistance is not None:
		barcodeDf['counts'] = barcodeDf.counts_x + barcodeDf.counts_y
		barcodeDf.drop(['counts_x', 'counts_y'], axis=1, inplace = True)
	else:	
		barcodeDf['reads'] = barcodeDf.reads_x + barcodeDf.reads_y
		barcodeDf['UMIs'] = barcodeDf.UMIs_x + barcodeDf.UMIs_y
		barcodeDf.drop(['reads_x', 'reads_y', 'UMIs_x', 'UMIs_y'], axis=1, inplace = True)
else:
	barcodeDf = barcodeDataFrames[0]


if args.starcodeDistance is not None:
	barcodeDftopN = barcodeDf.sort_values(by = 'counts', ascending = False).iloc[0:args.number]
else:
	barcodeDftopN = barcodeDf.sort_values(by = args.countMetric, ascending = False).iloc[0:args.number]

if args.mask is not None:
	#read mask sequence file. 
	with open(args.mask, 'r') as maskFile:
		header = maskFile.readline().strip('\n')
		maskSequence = maskFile.readline().strip('\n')
	barcodeDftopN['distanceToMask'] = barcodeDftopN.apply(lambda x: getMinHamming(x['barcode'], maskSequence), axis = 1)
	barcodeDfkeep = barcodeDftopN[barcodeDftopN['distanceToMask'] > args.numMask]
	barcodeDfdrop = barcodeDftopN[barcodeDftopN['distanceToMask'] <= args.numMask]
else:
	barcodeDfkeep = barcodeDftopN


barcodeDfkeep = barcodeDfkeep[barcodeDfkeep.barcode.str.len() >= args.minLength]
barcodeDfdrop = pd.concat([barcodeDfdrop, barcodeDfkeep[barcodeDfkeep.barcode.str.len() < args.minLength]], axis = 0)

#Reverse complement the barcode sequences and write to fasta files
barcodeDfkeep['reverse_complement'] = barcodeDfkeep.barcode.apply(reverse_complement)
barcodeDfkeep.reset_index(drop = True, inplace = True)
barcodeDfdrop['reverse_complement'] = barcodeDfdrop.barcode.apply(reverse_complement)
barcodeDfdrop.reset_index(drop = True, inplace = True)

if args.extend is not None:
	extension = args.extend
	barcodeDfkeep['reverse_complement'] = extension[0] + barcodeDfkeep['reverse_complement'] + extension[1]

if args.outFileName is not None:
	outFilePrefix = args.outFileName
else: 	
	outFilePrefix = samples[0]

if not os.path.exists("barcodeDesign/{}/fastaFiles/".format(outFilePrefix)):
	os.makedirs("barcodeDesign/{}/fastaFiles/".format(outFilePrefix))       

if args.outputMasked == True:
	if not os.path.exists("barcodeDesign/{}/fastaFiles/badMask/".format(outFilePrefix)):
		os.makedirs("barcodeDesign/{}/fastaFiles/badMask/".format(outFilePrefix)) 
	barcodeDfkeep.apply(lambda x: writeToFasta("barcodeDesign/{}/fastaFiles/".format(outFilePrefix), x.name, x.reverse_complement), axis = 1)
	barcodeDfdrop.apply(lambda x: writeToFasta("barcodeDesign/{}/fastaFiles/badMask/".format(outFilePrefix), x.name, x.reverse_complement), axis = 1)
else:
	barcodeDfkeep.apply(lambda x: writeToFasta("barcodeDesign/{}/fastaFiles/".format(outFilePrefix), x.name, x.reverse_complement), axis = 1)
















