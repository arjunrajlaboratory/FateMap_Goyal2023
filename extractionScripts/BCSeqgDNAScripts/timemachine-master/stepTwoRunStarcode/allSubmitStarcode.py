from argparse import ArgumentParser
import os, sys
import glob
import subprocess

parser = ArgumentParser()
parser.add_argument("experiment", help = "Specify the path to the experiment directory")
parser.add_argument("-d", "--distance", help = "Specify the Levenshtein distance for clustering/merging sequences", type = str, default = 8) 
parser.add_argument("--check_vector",  help = "Were full length barcodes extracted based on the vector sequence before and after the barcode, or only before?", default = "both", choices = ["both", "before"])
parser.add_argument("-c", "--countType", help = "Specify whether to use read counts or UMI counts", default = "Reads", choices = ["Reads", "UMIs"])
parser.add_argument("-r", "--clusterRatio", help = "Specify the minimum count ratio to cluster sequences within specified distance. Default is 5", type = str, default = 5)    
args = parser.parse_args()

starcodeScriptPath = "/project/arjunrajlab/timeMachine/repo/barcodeAnalysisScripts/stepTwoRunStarcode/submitStarcode.py"

samples = [os.path.basename(i) for i in glob.glob("{}/analyzed/*".format(args.experiment))]

#Run starcode on each sample file
for index, i in enumerate(samples):
	command = ["bsub", "-M", "32000", "-J", "starcode{}_d{}".format(index + 1, args.distance), "-o", "{}.starcode_d{}.stdout".format(i, args.distance), "-e", "{}.starcode_d{}.stderr".format(i, args.distance), \
				"python", starcodeScriptPath, args.experiment, i, "-d", args.distance, "--check_vector", args.check_vector, "--countType", args.countType, "--clusterRatio", args.clusterRatio]
	subprocess.call(command)
	#print command