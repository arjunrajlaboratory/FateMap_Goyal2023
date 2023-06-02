from argparse import ArgumentParser
import os, sys
import glob
import subprocess


#Command line parser
parser = ArgumentParser()
parser.add_argument("experiment", help = "Specify the path to the experiment directory")
parser.add_argument("sampleName", help = "Specify the name of sample directory containing the fastq.gz files")
parser.add_argument("-d", "--distance", help = "Specify the Levenshtein distance for clustering/merging sequences", type = str, default = 8) 
parser.add_argument("--check_vector",  help = "Were full length barcodes extracted based on the vector sequence before and after the barcode, or only before?", default = "both", choices = ["both", "before"])
parser.add_argument("-c", "--countType", help = "Specify whether to use read counts or UMI counts", default = "Reads", choices = ["Reads", "UMIs"])  
parser.add_argument("-r", "--clusterRatio", help = "Specify the minimum count ratio to cluster sequences within specified distance. Default is 5", type = str, default = 5) 
args = parser.parse_args()

#Create file extension
if args.countType == "Reads":
	countFileExtension = "readCountsOnly"
elif args.countType == "UMIs":
	countFileExtension = "UMICountsOnly"

if args.check_vector == "before":
	countFileExtension += "_liberal"

#First need to unzip barcode count file
samplePath = os.path.join(args.experiment, "analyzed", args.sampleName, 'extractedBarcodeData/{}_{}.gz'.format(args.sampleName, countFileExtension))
if os.path.isfile(samplePath):
	unzipCommand = ['gzip', '-d', samplePath]
	subprocess.call(unzipCommand)
else:
	samplePath = os.path.join(args.experiment, "analyzed", args.sampleName, 'extractedBarcodeData/{}_{}'.format(args.sampleName, countFileExtension))
	if os.path.isfile(samplePath) == False:
		print "Count file does not exist for sample {}".format(args.sampleName)
		sys.exit()

outFileDirectory = os.path.join(args.experiment, "analyzed", args.sampleName, 'starcode')		
if not os.path.exists(outFileDirectory):
        os.makedirs(outFileDirectory)

samplePath = os.path.join(args.experiment, "analyzed", args.sampleName, 'extractedBarcodeData/{}_{}'.format(args.sampleName, countFileExtension))
outfilePath = os.path.join(args.experiment, "analyzed", args.sampleName, 'starcode', '{}_clusteredBarcode{}_d{}.txt'.format(args.sampleName, args.countType, args.distance))
starcodeCommand = ['starcode', '-d', args.distance, '-t', str(4), '-r', args.clusterRatio, '-i', samplePath, '-o', outfilePath]

subprocess.call(starcodeCommand)