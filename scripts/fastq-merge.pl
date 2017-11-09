#!/usr/bin/env perl

use strict;
use warnings;

my $HELP = <<HELP;

  Find FASTQ files in a given directory (must have "_R1" or "_1" in file name).
  Extract sample names and paired reads based on file names (ignore lane numbers and barcodes).
  Combine multiple FASTQs so each sample ends up with only one R1 and R2 file.

  Usage:
    perl merge-fastqs.pl search_dir out_dir

  Arguments:
    search_dir   directory with original files
    out_dir      output directory for merged files

HELP

if (!$ARGV[1]) {
	die $HELP;
}

main();

# main subroutine
sub main {
	my $search_dir = $ARGV[0];
	my $out_dir = $ARGV[1];

	# convert dir from relative to absolute
	$search_dir = `readlink -f $search_dir`;
	chomp($search_dir);
	$out_dir = `readlink -f $out_dir`;
	chomp($out_dir);

	# check that given directories exist
	unless ( -d $search_dir ) {
		die "\n\n ERROR! $search_dir DOES NOT EXIST \n\n";
	}
	unless ( -d $out_dir ) {
		die "\n\n ERROR! $out_dir DOES NOT EXIST \n\n";
	}

	# find fastqs in given directory
	my $find_fastq_cmd_names = "-name '*_R1_0*.fastq.gz' -or -name '*_R1.fastq.gz' -or -name '*_1.fastq.gz'";
	my $find_fastq_cmd = "find -L $search_dir -maxdepth 3 -type f $find_fastq_cmd_names | LC_ALL=C sort";
	my @fastqs = `$find_fastq_cmd`;

	# process each fastq
	while (my $fastq_r1 = shift(@fastqs)) {
		chomp($fastq_r1);

		# check that R1 exists
		unless ( -e $fastq_r1 ) {
			die "\n\n ERROR! $fastq_r1 DOES NOT EXIST \n\n";
		}

		# generate R2 filename
		my $fastq_r2 = $fastq_r1;
		$fastq_r2 =~ s/(.*)_R1_0([0-9]+.fastq.gz)/${1}_R2_0${2}/;
		$fastq_r2 =~ s/(.*)_R1.fastq.gz/${1}_R2.fastq.gz/;
		$fastq_r2 =~ s/(.*)_1.fastq.gz/${1}_2.fastq.gz/;

		# blank if R2 does not exist
		unless ( -e $fastq_r2 ) {
			$fastq_r2 = "";
		}

		# blank if R2 is same as R1 (in case of not standard file name, for example)
		if ( $fastq_r1 eq $fastq_r2 ) {
			$fastq_r2 = "";
		}

		# extract sample name
		my $sample = $fastq_r1;
		# remove directory structure
		$sample =~ s/.*\///;
		# bcl2fastq2 format (with S sample number)
		$sample =~ s/_S[0-9]{1,3}_L00[0-9]_R1.*//;
		# bcl2fastq format with 2 barcodes
		$sample =~ s/_[ACTG]{6,}-[ACTG]{6,}_L00[0-9]_R1.*//;
		# bcl2fastq format with 1 barcode
		$sample =~ s/_[ACTG]{4,}_L00[0-9]_R1.*//;
		# no barcodes
		$sample =~ s/_L00[0-9]_R[12].*//;
		# no barcodes or lane
		$sample =~ s/_R[12].fastq.gz//;
		# no barcodes or lane
		$sample =~ s/_[12].fastq.gz//;

		# show found sample info
		print STDERR " SAMPLE : $sample \n";
		print STDERR "  FASTQ R1 : $fastq_r1 \n";
		print STDERR "  FASTQ R2 : $fastq_r2 \n";

		# file names after merging
		my $fastq_r1_merged = "${out_dir}/${sample}_R1.fastq.gz";
		my $fastq_r2_merged = "${out_dir}/${sample}_R2.fastq.gz";

		# merge without overwriting (">>" instead of ">")
		print STDERR "  MERGED R1 : $fastq_r1_merged \n";
		my $merge_cmd = "cat $fastq_r1 >> $fastq_r1_merged";
		print STDERR "  CMD : $merge_cmd \n";
		system($merge_cmd);

		# repeat merge for R2 if present
		if ( -e $fastq_r2 ) {
			print STDERR "  MERGED R2 : $fastq_r2_merged \n";
			$merge_cmd = "cat $fastq_r2 >> $fastq_r2_merged";
			print STDERR "  CMD : $merge_cmd \n";
			system($merge_cmd);
		}

		sleep(1);

	}

}



# end
