#!/usr/bin/perl

# run from within repo directory.
# takes relative paths from the repo directory as arguments.

use File::Spec;

print "experiment\tsampleID\tgene_id\tcounts\n";

while (@ARGV) {
	$fullPath = shift @ARGV;
	($volume,$directories,$file) = File::Spec->splitpath($fullPath);

	@dirs = File::Spec->splitdir($directories);

	$experiment = shift @dirs;
	shift @dirs;
	$sampleID = shift @dirs;

	open FILE, "< $fullPath";

	for (<FILE>) {
		print "$experiment\t$sampleID\t";
		print;
	}
}