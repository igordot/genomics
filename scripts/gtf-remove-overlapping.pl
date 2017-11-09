#!/usr/bin/perl

use strict;

my $HELP = <<HELP;

  REMOVE GENES FROM GTF FILE WITH OVERLAPPING REGIONS

  Usage: gtf-remove-overlapping.pl in.gtf
         (will create a new filtered file in the same dir)
HELP

if(!$ARGV[0]){
	die $HELP;
}

main();

# main subroutine
sub main {
	my $gtf = $ARGV[0];

	# check that the gtf file exists
	unless ( -s $gtf ) {
		print "\n\n ERROR! $gtf DOES NOT EXIST \n\n";
		die;
	}

	# get file stats
	my $gtf_lines = `cat $gtf | cut -f 9 | wc -l`;
	my $gtf_transcripts = `cat $gtf | cut -f 9 | sort | uniq | wc -l`;
	print " original gtf lines: $gtf_lines";
	print " original gtf transcripts: $gtf_transcripts";

	# keep only cds entries and sort gtf
	my $gtf_sorted = "${gtf}.sort.tmp";
	my $sort_cmd = "grep 'CDS' $gtf | sort -k1,1 -k4,4n -k5,5n -k9,9n > $gtf_sorted";
	system($sort_cmd);
	sleep(1);

	# scan the gtf file one line at a time checking for overlapping transcripts
	my ($prev_chr, $prev_pos0, $prev_pos1, $prev_transcript);
	my @bad_transcripts;
	open(SORTED, "<", $gtf_sorted);
	while (<SORTED>) {
		chomp;
		my ($chr, $src, $feat, $pos0, $pos1, $score, $strand, $frame, $attr) = split(/\t/);

		# extract just the transcript id from attributes column
		my $transcript = $attr;
		$transcript =~ s/.*transcript_id\s"(.+?)".*/$1/;

		# check if current region is same chr as previous region but transcript is different
		if ( ( $chr eq $prev_chr ) && ( $transcript ne $prev_transcript ) ) {
			# check if previous end pos is larger than current start pos
			if ( $pos0 < $prev_pos1 ) {
				# flag previous transcript if not flagged already
				unless ($prev_transcript ~~ @bad_transcripts) {
					push (@bad_transcripts, $prev_transcript);
					#print "$prev_transcript \n";
				}
				# flag current transcript if not flagged already
					unless ($transcript ~~ @bad_transcripts) {
					push (@bad_transcripts, $transcript);
					#print "$transcript \n";
				}
			}
		}

		#print " $i $gtf_array[$i][0] $gtf_array[$i][3] $gtf_array[$i][4] \n";

		# update info of previous entry
		$prev_chr = $chr;
		$prev_pos0 = $pos0;
		$prev_pos1 = $pos1;
		$prev_transcript = $transcript;


		#for my $row (@gtf_array) {
		#	print "@$row[0]\t@$row[1]\t@$row[2]\n";
		#}
	}
	close(SORTED);

	# delete sorted cds-only gtf created at the beginning
	system("rm -f $gtf_sorted");

	# count number of overlapping transcripts
	my $bad_transcript_count = scalar(@bad_transcripts);
	print " overlapping transcripts: $bad_transcript_count \n";

	# create new gtf for overlapping and non-overlapping transcripts
	open(UNIQUE, ">", "${gtf}.unique.gtf");
	open(OVERLAPPING, ">", "${gtf}.overlapping.gtf");

	# process original gtf and sort each element based on transcript being one of overlapping transcripts
	open(GTF, "<", $gtf);
	while (<GTF>) {
		chomp;
		my ($chr, $src, $feat, $pos0, $pos1, $score, $strand, $frame, $attr) = split(/\t/);

		# extract just the transcript id from attributes column
		my $transcript = $attr;
		$transcript =~ s/.*transcript_id\s"(.+?)".*/$1/;

		# check if transcript is one of overlapping transcripts
		if ($transcript ~~ @bad_transcripts) {
			print OVERLAPPING "$_\n";
		}
		else {
			print UNIQUE "$_\n";
		}
	}
	close(GTF);

	close(UNIQUE);
	close(OVERLAPPING);

}



