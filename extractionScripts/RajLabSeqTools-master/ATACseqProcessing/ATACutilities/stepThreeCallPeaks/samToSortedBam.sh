#!
# run from within "repo" directory

EXPERIMENT=$1
SAMPLEID=$2

ALIGNMENT_TOOL_NAME=bowtie2

if [[ ! -z "$3" ]]; then
    ALIGNMENT_TOOL_NAME="$3"
fi

COMMANDNAME=samtools

if [ ! -d $EXPERIMENT/analyzed/$SAMPLEID/log ]; then
    mkdir $EXPERIMENT/analyzed/$SAMPLEID/log
fi

JOURNAL=$EXPERIMENT/analyzed/$SAMPLEID/log/$(date +%Y-%m-%d_%H-%M).$COMMANDNAME.log

BEDTOOLS_LOC="/Users/sydneyshaffer/ATACseqTools/bedtools2-master/bin"

chromSizes="/Users/sydneyshaffer/ATACseqTools/hg19_chromosomeSizes.txt"

# opens SAM file, sorts, coverts to BAM file
COMMAND="samtools view -bS \
    $EXPERIMENT/analyzed/$SAMPLEID/$ALIGNMENT_TOOL_NAME/$SAMPLEID.sam \
    | samtools sort - \
    $EXPERIMENT/analyzed/$SAMPLEID/$ALIGNMENT_TOOL_NAME/$SAMPLEID.sorted"

# writes only the reads that mapped to file
COMMAND_MAPPED="samtools view -b -F 4 \
	$EXPERIMENT/analyzed/$SAMPLEID/$ALIGNMENT_TOOL_NAME/$SAMPLEID.sorted.bam > \
	$EXPERIMENT/analyzed/$SAMPLEID/$ALIGNMENT_TOOL_NAME/$SAMPLEID.sorted.mapped.bam"

# removed duplicate reads
COMMAND_DUPS="samtools rmdup $EXPERIMENT/analyzed/$SAMPLEID/$ALIGNMENT_TOOL_NAME/$SAMPLEID.sorted.mapped.bam \
	$EXPERIMENT/analyzed/$SAMPLEID/$ALIGNMENT_TOOL_NAME/$SAMPLEID.sorted.mapped.rmDups.bam"

# make index file
INDEXCOMMAND="samtools index \
    $EXPERIMENT/analyzed/$SAMPLEID/$ALIGNMENT_TOOL_NAME/$SAMPLEID.sorted.mapped.rmDups.bam"
    
# Remove chromsome Y and mito reads
COMMAND_RM="samtools view -b \
	$EXPERIMENT/analyzed/$SAMPLEID/$ALIGNMENT_TOOL_NAME/$SAMPLEID.sorted.mapped.rmDups.bam \
	chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 \
	chr17 chr18 chr19 chr20 chr21 chr22 chrX >\
	$EXPERIMENT/analyzed/$SAMPLEID/$ALIGNMENT_TOOL_NAME/$SAMPLEID.sorted.mapped.rmDups.rmMandY.bam"
	
# Make bedgraph file
BEDGRAPH="$BEDTOOLS_LOC/genomeCoverageBed -bg -ibam \
	$EXPERIMENT/analyzed/$SAMPLEID/$ALIGNMENT_TOOL_NAME/$SAMPLEID.sorted.mapped.rmDups.rmMandY.bam \
	-split -g $chromSizes >\
	$EXPERIMENT/analyzed/$SAMPLEID/$ALIGNMENT_TOOL_NAME/$SAMPLEID.bedGraph"

echo "Starting..." >> $JOURNAL

date >> $JOURNAL
echo $COMMAND >> $JOURNAL
eval $COMMAND

date >> $JOURNAL
echo $COMMAND_MAPPED >> $JOURNAL
eval $COMMAND_MAPPED

date >> $JOURNAL
echo $COMMAND_DUPS >> $JOURNAL
eval $COMMAND_DUPS

date >> $JOURNAL
echo $INDEXCOMMAND >> $JOURNAL
eval $INDEXCOMMAND

date >> $JOURNAL
echo $COMMAND_RM >> $JOURNAL
eval $COMMAND_RM

date >> $JOURNAL
echo $BEDGRAPH >> $JOURNAL
eval $BEDGRAPH

date >> $JOURNAL
echo "Done" >> $JOURNAL
