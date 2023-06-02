from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip
import regex as re
from collections import Counter
import os, glob
import numpy as np

def parseBarcodeFastQ_both(inFileNames, staggerLength, minPhred):
	"""Function to parse FASTQ file and create dictionary of key:value pairs. Will only extract reads that contain the vector sequence before and after the barcode, allowing up to 4 or 5 mismatches respectively. 
	Each key is a unique barcode sequence. Each value is a list containing the Phredscore associated with the barcode
	and the UMIs (first 4-6 bases) associated with the read. Function will also output the list of reads that are missing vector sequence before the barcode, missing the vector sequence after the barcode, have a bad Qscore before the barcode or contain a short UMI or primer binding sequence. 
	"""
	barcode_dict = {} 
	missingBeforeBarcode = []
	missingAfterBarcode = []
	badQscore = []
	badLength = []
	vectorBeforeBarcode = re.compile(r'(?e)(TCGACTAAACGCGCTACTTGAT){e<=4}') #vector primer site as reg expression. Allow up to 4 mismatches.
	vectorAfterBarcode = re.compile(r'(?e)(?r)(ATCCTACTTGTACAGCTCGT){e<=5}') #vector sequence after barcode as reg expression. Allow up to 5 mismatches and search from end of string first.
	for fastQ_file in inFileNames:
	#Loop over all fastQ files/lanes associated with each sample
		with gzip.open(fastQ_file, 'rb') as fastq:
			for seq_record in FastqGeneralIterator(fastq): #Uses BioPythons FastqGeneralIterator to parse each read.
				if vectorBeforeBarcode.search(seq_record[1]) and vectorAfterBarcode.search(seq_record[1][80:]): #Check that the read contains the vector surrounding the barcode.
					spanBeforeBarcode = vectorBeforeBarcode.search(seq_record[1]).span()
					spanAfterBarcode = vectorAfterBarcode.search(seq_record[1]).span() 
					if sum([i < minPhred for i in seq_record[2][spanBeforeBarcode[0]:spanBeforeBarcode[1]]]) >= 5: #Skip reads where >=5 positions in GFP primer site have phredscore < 14.   
						badQscore.append(seq_record)
					elif (spanBeforeBarcode[1] - staggerLength) > 30 or (spanBeforeBarcode[0]-staggerLength) < 4: #Skip positions with UMI shorter than 4 bases or UMI+GFP longer than 30 bases. 
						badLength.append(seq_record)
					else:
						barcode = seq_record[1][spanBeforeBarcode[1]:spanAfterBarcode[0]]
						if barcode not in barcode_dict:
							barcode_dict[barcode] = [seq_record[2][spanBeforeBarcode[1]:spanAfterBarcode[0]]]
							barcode_dict[barcode].append(seq_record[1][0:spanBeforeBarcode[0]-staggerLength])
						elif barcode in barcode_dict:
							barcode_dict[barcode].append(seq_record[1][0:spanBeforeBarcode[0]-staggerLength])
				elif vectorBeforeBarcode.search(seq_record[1]) and not vectorAfterBarcode.search(seq_record[1][80:]):
					missingAfterBarcode.append(seq_record)
				elif vectorAfterBarcode.search(seq_record[1][80:]) and not vectorBeforeBarcode.search(seq_record[1]):
					missingBeforeBarcode.append(seq_record)
		print "Completed file: " + fastQ_file
	return barcode_dict, missingBeforeBarcode, missingAfterBarcode, badQscore, badLength

def parseBarcodeFastQ_before(inFileNames, staggerLength, barcodeLength, minPhred):
	"""Function to parse FASTQ file and create dictionary of key:value pairs. Will only extract reads that contain the vector sequence before the barcode allowing up to 4 mismatches. 
	Each key is a unique barcode sequence. Each value is a list containing the Phredscore associated with the barcode
	and the UMIs (first 4-6 bases) associated with the read. Function will also output the list of reads that are missing vector sequence before the barcode, missing the vector sequence after the barcode, have a bad Qscore before the barcode or contain a short UMI or primer binding sequence.
    """
	barcode_dict = {} 
	missingBeforeBarcode = []
	badQscore = []
	badLength = []
	vectorBeforeBarcode = re.compile(r'(?e)(TCGACTAAACGCGCTACTTGAT){e<=4}') #vector primer site as reg expression. Allow up to 4 mismatches.
	vectorAfterBarcode = re.compile(r'(?e)(?r)(ATCCTACTTGTACAGCTCGT){e<=5}') #vector sequence after barcode as reg expression. Allow up to 5 mismatches and search from end of string first.
	for fastQ_file in inFileNames:
	#Loop over all fastQ files/lanes associated with each sample
		with gzip.open(fastQ_file, 'rb') as fastq:
			for seq_record in FastqGeneralIterator(fastq): #Uses BioPythons FastqGeneralIterator to parse each read.
				if vectorBeforeBarcode.search(seq_record[1]): #Check that the read contains the vector sequence before the barcode.
					spanBeforeBarcode = vectorBeforeBarcode.search(seq_record[1]).span()
					if sum([i < minPhred for i in seq_record[2][spanBeforeBarcode[0]:spanBeforeBarcode[1]]]) >= 5: #Skip reads where >=5 positions in GFP primer site have phredscore < 14.   
						badQscore.append(seq_record)
					elif (spanBeforeBarcode[1] - staggerLength) > 30 or (spanBeforeBarcode[0]-staggerLength) < 4: #Skip positions with UMI shorter than 4 bases or UMI+GFP longer than 30 bases. 
						badLength.append(seq_record)
					else:
						barcode = seq_record[1][spanBeforeBarcode[1]:spanBeforeBarcode[1] + barcodeLength]
						if barcode not in barcode_dict:
							barcode_dict[barcode] = [seq_record[2][spanBeforeBarcode[1]:spanBeforeBarcode[1] + barcodeLength]]
							barcode_dict[barcode].append(seq_record[1][0:spanBeforeBarcode[0]-staggerLength])
						elif barcode in barcode_dict:
							barcode_dict[barcode].append(seq_record[1][0:spanBeforeBarcode[0]-staggerLength])
				elif vectorAfterBarcode.search(seq_record[1][80:]) and not vectorBeforeBarcode.search(seq_record[1]):
					missingBeforeBarcode.append(seq_record)
		print "Completed file: " + fastQ_file
	return barcode_dict, missingBeforeBarcode, badQscore, badLength

def writeOutFileUMIs(barcode_dict, outFileName):
    """Function to write barcode dictionary to gzipped tab delimited file. Each line contains
    a unique barcode, it's phredscore, and all associated UMIs. 
    """
    with gzip.open(outFileName, 'wb') as out_file:
    	for barcode in barcode_dict:
            out_file.write(barcode)
            out_file.write("\t" + "\t".join(barcode_dict[barcode]))
            out_file.write("\n")

def writeOutFileBadSeqRecord(badSeqList, outFileName):
    """Function to write a list of bad sequence records to gzipped tab delimited file. 
    """
    with gzip.open(outFileName, 'wb') as out_file:
    	for seqRecord in badSeqList:
            out_file.write("\t".join(seqRecord) + "\n")

def countUMIs(barcode_dictionary):
    """Function to count the number of reads and the number of UMIs associated with each barcode. 
    Output is dictionary with the latter as the key and the former as the associated values
    """
    new_dict = {}
    for i in barcode_dictionary:
    	new_dict[i] = (sum(Counter(barcode_dictionary[i][1:]).values()), len(Counter(barcode_dictionary[i][1:]).keys()))
    return new_dict

def writeOutFileBarcodeCounts(barcode_dict_summary, outFileName):
    """Function to write reads counts and UMIs associated with each barcode
    to gzipped, tab-delimited file.
    """
    with gzip.open(outFileName, 'wb') as out_file:
        for barcode in barcode_dict_summary:
            out_file.write(barcode)
            out_file.write("\t" + "\t".join(map(str,barcode_dict_summary[barcode])))
            out_file.write("\n")    

def writeOutFileBarcodeReadCounts(barcode_dict_summary, outFileName):
    """Function to write reads counts and UMIs associated with each barcode
    to gzipped, tab-delimited file.
    """
    with gzip.open(outFileName, 'wb') as out_file:
        for barcode in barcode_dict_summary:
            out_file.write(barcode)
            out_file.write("\t" + str(barcode_dict_summary[barcode][0]))
            out_file.write("\n")    

def writeOutFileBarcodeUMICounts(barcode_dict_summary, outFileName):
    """Function to write UMIs counts associated with each barcode
    to gzipped, tab-delimited file.
    """
    with gzip.open(outFileName, 'wb') as out_file:
        for barcode in barcode_dict_summary:
            out_file.write(barcode)
            out_file.write("\t" + str(barcode_dict_summary[barcode][1]))
            out_file.write("\n")    

