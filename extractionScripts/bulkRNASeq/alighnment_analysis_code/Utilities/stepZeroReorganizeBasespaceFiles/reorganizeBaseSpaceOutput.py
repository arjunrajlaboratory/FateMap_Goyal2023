import os, glob, shutil, re
from argparse import ArgumentParser

#Command line parser
parser = ArgumentParser()
parser.add_argument("path", help = "Specify the path containing the sample files.", type = str)
parser.add_argument("-n", "--names", help = "Option to specify one or more sample names. If specified, will use this name find directories containing FASTQ files from each lane and move these files to new directory 'name'", nargs = "+", type = list)
args = parser.parse_args()

os.chdir(args.path)

oldSampleDirectories = glob.glob("*L00*")
if args.names is None:
	laneRegex = re.compile('L\d+')
	sampleNames = list(set([i[0:laneRegex.search(i).span()[0]-1] for i in oldSampleDirectories]))
else:
	sampleNames = args.names

for sample in sampleNames:
	if not os.path.isdir(sample):
		os.makedirs(sample)
	sampleLanes = glob.glob("{}*_L00*".format(sample))
	for i in sampleLanes:
		for file in glob.glob("{}/*.fastq.gz".format(i)):
			shutil.move(file, "{}/".format(sample))
		shutil.rmtree(i)	


