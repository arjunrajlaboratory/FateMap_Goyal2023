import os, glob, shutil, re
from argparse import ArgumentParser

#Command line parser
parser = ArgumentParser()
parser.add_argument("path", help = "Specify the path containing the sample files.", type = str)
args = parser.parse_args()

os.chdir(args.path)

oldSampleDirectories = glob.glob("*L00*")
laneRegex = re.compile('L\d+')
sampleNames = list(set([i[0:laneRegex.search(i).span()[0]-1] for i in oldSampleDirectories]))

for sample in sampleNames:
	if not os.path.isdir(sample):
		os.makedirs(sample)
	sampleLanes = glob.glob("{}_*L00*".format(sample))
	for i in sampleLanes:
		for file in glob.glob("{}/*.fastq.gz".format(i)):
			shutil.move(file, "{}/".format(sample))
		shutil.rmtree(i)	


