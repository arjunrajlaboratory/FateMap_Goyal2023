#!/bin/bash
# run this script before running any submission scripts!
# this script allows you to use samtools, ngstools, etc.

## Do not change this section of code unless you want to change software versions!
# within an interactive node (type "bsub -Is bash" to get to one), 
# you can see the list of available software versions by typing "module avail"
module load FastQC-0.11.2
module load samtools-0.1.19
module load ngsutils-0.5.7
module load STAR/2.5.2a
module load python/2.7.9                               # Loading python 2 is necessary because attempting to run with python 3 will lead to import errors at the HTSeq step.
STARFLAGS="--readFilesCommand zcat"
genomeDirSTAR="/home/apps/STAR/indexes/hg19"    # This file contains an index used by the STAR aligner. Change the index if you're not using the hg38 reference genome. For hg19, use /home/apps/STAR/indexes/hg19
gtfFile="/project/arjunrajlab/resources/htseq/hg19/hg19.gtf"      # This file contains transcript information for hg38 genes. Change it to a different file for a different reference genome. For hg19, use /project/arjunrajlab/resources/htseq/hg19/hg19.gtf
#######################

## Update these variables below to your project name, experiment name, number of samples,
#  code home directory (where your "rajlabseqtools/Utilities" folder is. if it's in your
#  home directory (e.g. /home/esanford), you can leave the "~" symbol)
PROJECT="mapseq"
EXPERIMENT="BulkColonySeq_20200220"
RAWDATA_DIRECTORY="/project/arjunrajlab/mapseq/repo/BulkColonySeq_20200220/raw_fastqz/run3"
N_SAMPLES=67
PAIRED_OR_SINGLE_END_FRAGMENTS="paired"  # this variable must be "single" or "paired". change to "paired" if your reads... are paired.
codeHomeDir=~                            # "~" is a shortcut for your home directory. alternatively you can use /home/<YOUR_PMACS_USERNAME>



