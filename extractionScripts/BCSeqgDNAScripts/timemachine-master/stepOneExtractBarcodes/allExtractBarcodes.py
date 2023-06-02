import os, subprocess
import csv
from argparse import ArgumentParser
import glob

#Command line parser
parser = ArgumentParser()
parser.add_argument("experiment", help = "Specify the path to the experiment directory")
parser.add_argument("--staggerFile", help = "Specify the name of the file containing information on length of stagger for each sample")
parser.add_argument("-r", "--includeReads", help = "If specified, output additional tables with barcodes and only read counts or UMI counts. Otherwise outputs only one table with both counts", action = 'store_true')
parser.add_argument("--check_vector", help = "Option to check vector sequence before or on both sides of the barcode sequence.", default = "both", choices = ["both", "before"]) 
parser.add_argument("-l", "--length", help = "If check_vector before specified, input here your desired barcode length. Default is 100", default = "100", type = str) 
parser.add_argument("-Q", "--minPhred", help = "Specify the minimum phredscore required to include a readout. Filters reads with more than 5 bases before the barcode with low phredscore.", default = "14", type = str) 
parser.add_argument("-e", "--excludedReads", help = "If specified, output txt.gz files containing reads excluded from the UMI and count files.", action = 'store_true')
args = parser.parse_args()

extractionScriptPath = "/project/arjunrajlab/timeMachine/repo/barcodeAnalysisScripts/stepOneExtractBarcodes/timeMachineParseFastQ150SEreads_v4.py"
#loadPythonEnvPath = '/barcodeAnalysisScripts/stepOneExtractBarcodes/set_python_pmacs_env'
#additionalArguments = ["--includeReads", "--check_vector", "both", "--excludedReads"]

#Move to experiment directory
os.chdir(args.experiment)
if args.staggerFile is not None:
	samples = []
	staggers = []
	with open(args.staggerFile, 'r') as file:
		tmp = csv.reader(file)
		for line in tmp:
			samples.append(line[0])
			staggers.append(str(line[1]))
else:
	samples = [os.path.basename(i) for i in glob.glob("raw/*")]
	staggers = ['0'] * len(samples)

#Format additional arguments
if args.includeReads:
	additionalArguments = ["--includeReads"]
else:
	additionalArguments = []

additionalArguments.extend(["--check_vector", args.check_vector])

if args.check_vector == "before":
	additionalArguments.extend(["--length", args.length])

additionalArguments.extend(["--minPhred", args.minPhred])

if args.excludedReads:
	additionalArguments.append("--excludedReads")

print additionalArguments
#Extract barcode for each sample file
for index, i in enumerate(samples):
	command = ["bsub",  "-J", "exractBarcodes{}".format(index + 1), "-o", "{}.extractBarcodes.stdout".format(i), "-e", "{}.extractBarcodes.stderr".format(i), \
				"python", extractionScriptPath, i, "-s", staggers[index]] \
				+ additionalArguments
	print command
	subprocess.call(command)
	#print command