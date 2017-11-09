#!/usr/bin/env perl

use strict;
use warnings;
use v5.10;
use File::Basename;
use Bio::SeqIO;

my $HELP = <<HELP;

  Predict MHC-I binding for ANNOVAR-annotated point mutations using IEDB MHC class I binding prediction tool.
  Protein coding changes are annotated using ANNOVAR.
  Binding predictions are calculated by IEDB MHC class I binding prediction tool using consensus method.

  Usage:
    perl mut-mhc-i-binding.pl base_name avinput_file annovar_scripts_dir annovar_refgene_txt iedb_mhc_i_dir

  Arguments:
    base_name             output prefix
    avinput_file          ANNOVAR avinput file - 5 required cols (chr, start, end, ref, alt) and then optional cols
    annovar_scripts_dir   ANNOVAR scripts directory
    annovar_refgene_txt   ANNOVAR refGene.txt reference file
    iedb_mhc_i_dir        IEDB MHC_I binding tool directory (tested with version 2.17)

HELP

if (!$ARGV[4]) {
	die $HELP;
}

main();

# main subroutine
sub main {
	my $base_name = $ARGV[0];
	my $avinput = $ARGV[1];
	my $annovar_scripts_dir = $ARGV[2];
	my $annovar_refgene_txt = $ARGV[3];
	my $iedb_mhc_i_dir = $ARGV[4];

	# check that inputs exist
	unless ( -e $avinput ) {
		die "\n\n ERROR: $avinput does not exist \n\n";
	}
	unless ( -d $annovar_scripts_dir ) {
		die "\n\n ERROR: $annovar_scripts_dir does not exist \n\n";
	}
	unless ( -e $annovar_refgene_txt ) {
		die "\n\n ERROR: $annovar_refgene_txt does not exist \n\n";
	}
	unless ( -d $iedb_mhc_i_dir ) {
		die "\n\n ERROR: $iedb_mhc_i_dir does not exist \n\n";
	}

	# check that the ANNOVAR and IEDB scripts exist
	my $annovar_annotate_pl = "${annovar_scripts_dir}/annotate_variation.pl";
	unless ( -e $annovar_annotate_pl ) {
		die "\n\n ERROR: $annovar_annotate_pl does not exist \n\n";
	}
	my $annovar_coding_change_pl = "${annovar_scripts_dir}/coding_change.pl";
	unless ( -e $annovar_coding_change_pl ) {
		die "\n\n ERROR: $annovar_coding_change_pl does not exist \n\n";
	}
	my $predict_binding_py = "${iedb_mhc_i_dir}/src/predict_binding.py";
	unless ( -e $predict_binding_py ) {
		die "\n\n ERROR: $predict_binding_py does not exist \n\n";
	}

	# peptide length for binding predictions (and preparing binding prediction input sequences)
	my $peptide_length = 9;

	say "annotate avinput";
	my $evf = annotate_avinput($base_name, $avinput, $annovar_annotate_pl, $annovar_refgene_txt);

	say "get coding change";
	my $cod_change_fa = get_coding_change($base_name, $evf, $annovar_coding_change_pl, $annovar_refgene_txt);

	say "format coding change";
	my $mutpad_fa = format_coding_change_fa($base_name, $cod_change_fa, $peptide_length);

	say "predict binding";
	my $bindpred_txt = predict_mhc_i_binding($base_name, $mutpad_fa, $peptide_length, $predict_binding_py);

	say "annotate binding predictions";
	my $annot_txt = annotate_binding_predictions($base_name, $bindpred_txt, $evf);
	say "created $annot_txt";

}

# annotate avinput
sub annotate_avinput {
	my $base_name = $_[0];
	my $avinput = $_[1];
	my $annovar_annotate_pl = $_[2];
	my $annovar_refgene_txt = $_[3];

	my $genome = basename($annovar_refgene_txt);
	$genome =~ s/_refGene.txt//;
	my $annovar_ref_dir = dirname($annovar_refgene_txt);

	# run ANNOVAR annotate_variation.pl
	my $cmd = "perl $annovar_annotate_pl";
	$cmd .= " --geneanno --dbtype refGene --buildver $genome --outfile $base_name";
	$cmd .= " $avinput $annovar_ref_dir";
	system $cmd;

	# can only specify output file prefix
	my $evf = "${base_name}.exonic_variant_function";

	# confirm that exonic_variant_function file generated
	unless ( -e $evf ) {
		die "\n\n ERROR: $evf DOES NOT EXIST \n\n";
	}

	# clean up
	unlink "${base_name}.log";
	unlink "${base_name}.variant_function";

	return $evf;
}

# generate coding change FASTA
sub get_coding_change {
	my $base_name = $_[0];
	my $evf = $_[1];
	my $annovar_coding_change_pl = $_[2];
	my $annovar_refgene_txt = $_[3];

	my $annovar_refgene_fa = $annovar_refgene_txt;
	$annovar_refgene_fa =~ s/_refGene.txt/_refGeneMrna.fa/;
	my $out_file = "${base_name}.codchange.fa";

	# run ANNOVAR coding_change.pl
	my $cmd = "perl $annovar_coding_change_pl";
	$cmd .= " --includesnp --onlyAltering";
	$cmd .= " $evf $annovar_refgene_txt $annovar_refgene_fa";
	$cmd .= " > $out_file";
	system $cmd;

	# confirm that coding change file generated
	unless ( -e $out_file ) {
		die "\n\n ERROR: $out_file DOES NOT EXIST \n\n";
	}

	return $out_file;
}

# format coding change FASTA for IEDB MHC-I binding predictions
sub format_coding_change_fa {
	my $base_name = $_[0];
	my $coding_change_fa = $_[1];
	my $peptide_length = $_[2];

	# sequence padding (peptide length)
	my $padding = $peptide_length - 1;

	my $out_file = "${base_name}.mutpad.${peptide_length}.fa";

	# process each record
	my $seqio_in = Bio::SeqIO->new(-file => $coding_change_fa, -format => 'fasta');
	my $seqio_out = Bio::SeqIO->new(-file => ">$out_file", -format => 'fasta');
	my $seq_num = 0;
	while ( my $seq_in = $seqio_in->next_seq() ) {
		# extract relevant parts
		my $seq_id = $seq_in->id;
		my $seq_desc = $seq_in->desc;
		$seq_desc =~ s/\s+$//;
		my $sequence = $seq_in->seq;

		# process only altered sequences and single amino acid substitutions
		if ($seq_desc !~ m/silent/ && $seq_desc =~ m/\sp\.\w\d+\w\s/) {
			# extract mutation details and sequence around the altered amino acid
			my ($tx_name, $aa_num, $aa_from, $aa_to, $seq_padded) = ("ERR", "ERR", "ERR", "ERR", "ERR");
			if ($seq_desc =~ m/(.+?)\s.+?\(position\s(\d+)\schanged\sfrom\s(\w)\sto\s(\w)\)/) {
				$seq_num++;
				$tx_name = $1;
				$aa_num = $2;
				$aa_from = $3;
				$aa_to = $4;
				# adjust offsets if amino acid is close to the beginning of sequence
				my $substr_start = $aa_num - 1 - $padding;
				$substr_start = $substr_start < 0 ? 0 : $substr_start;
				my $substr_length = $aa_num > $padding ? $padding + $padding + 1 : $aa_num + $padding;
				# get sequence (0-based)
				$seq_padded = substr($sequence, $substr_start, $substr_length);
				# remove trailing asterisk (stop codon)
				$seq_padded =~ s/\*$//;
				# say "${seq_id}\t${tx_name}\t${aa_from}${aa_num}${aa_to}\t${seq_padded}";

				# create output FASTA record and write it (binding predictions will just assign consecutive numbers)
				my $new_id = "${seq_num}|${seq_id}|${tx_name}|${aa_from}${aa_num}${aa_to}|${seq_padded}";
				my $seq_out = Bio::Seq->new(-seq => $seq_padded, -id => $new_id);
				$seqio_out->write_seq($seq_out);
			}
			else {
				# say "skip";
			}
		}
	}

	# confirm that padded mutations FASTA file generated
	unless ( -e $out_file ) {
		die "\n\n ERROR: $out_file DOES NOT EXIST \n\n";
	}

	return $out_file;
}

# run IEDB MHC-I Binding Predictions
sub predict_mhc_i_binding {
	my $base_name = $_[0];
	my $mutpad_fa = $_[1];
	my $peptide_length = $_[2];
	my $predict_binding_py = $_[3];

	my $raw_out_file = "${base_name}.bindpred.${peptide_length}.raw.txt";
	my $out_file = "${base_name}.bindpred.${peptide_length}.txt";

	# delete output if already exists
	if ( -e $raw_out_file ) {
		unlink $raw_out_file;
	}

	# my @alleles = ('HLA-A*01:01');
	# my @alleles = ('HLA-A*01:01', 'HLA-A*02:01');

	# a reference panel of 27 alleles
	# http://help.iedb.org/hc/en-us/articles/114094151851
	my @alleles = ('HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*02:03', 'HLA-A*02:06', 'HLA-A*03:01', 'HLA-A*11:01',
		'HLA-A*23:01', 'HLA-A*24:02', 'HLA-A*26:01', 'HLA-A*30:01', 'HLA-A*30:02', 'HLA-A*31:01','HLA-A*32:01',
		'HLA-A*33:01', 'HLA-A*68:01', 'HLA-A*68:02', 'HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*15:01', 'HLA-B*35:01',
		'HLA-B*40:01', 'HLA-B*44:02', 'HLA-B*44:03', 'HLA-B*51:01', 'HLA-B*53:01', 'HLA-B*57:01', 'HLA-B*58:01');

	# make predictions for each allele
	foreach (@alleles) {
		say "running predictions for allele $_";

		# predict_binding.py command
		# ./src/predict_binding [method] [mhc] [peptide_length] [input_file]
		my $predict_binding_cmd = "$predict_binding_py consensus $_ $peptide_length $mutpad_fa >> $raw_out_file";
		system $predict_binding_cmd;
	}

	# confirm that binding predictions file generated
	unless ( -e $raw_out_file ) {
		die "\n\n ERROR: $raw_out_file DOES NOT EXIST \n\n";
	}

	# header for the combined file
	my $bindpred_header = `cat $raw_out_file | head -1 | cut -f 1,5-`;
	system "printf \"line_id\ttranscript_id\tp_change\taa_padded\t${bindpred_header}\" > $out_file";

	# clean up input files for joining
	system "cat $mutpad_fa | grep '^>' | cut -c 2- | tr '|' '\t' | LC_ALL=C sort -k1,1 > ${mutpad_fa}.txt";
	system "cat $raw_out_file | grep -v '^allele' | cut -f 1,2,5- | LC_ALL=C sort -k2,2 > ${raw_out_file}.tmp";

	# join, remove seq col
	my $join_cmd = 'LC_ALL=C join -t $\'\t\' -a 1 -a 2 -1 1 -2 2';
	$join_cmd .= " ${mutpad_fa}.txt";
	$join_cmd .= " ${raw_out_file}.tmp";
	$join_cmd .= " | cut -f 2-";
	$join_cmd .= " >> $out_file";
	system $join_cmd;

	# confirm that binding predictions file generated
	unless ( -e $out_file ) {
		die "\n\n ERROR: $out_file DOES NOT EXIST \n\n";
	}

	# clean up
	unlink $mutpad_fa;
	unlink "${mutpad_fa}.txt";
	unlink "${raw_out_file}.tmp";

	return $out_file;
}

# combine binding predictions table with variant annotations
sub annotate_binding_predictions {
	my $base_name = $_[0];
	my $bindpred_txt = $_[1];
	my $evf = $_[2];

	my $out_file = "${base_name}.bindpred.annot.txt";

	# header for the combined file
	my $bindpred_header = `cat $bindpred_txt | head -1 | cut -f 2-`;
	system "printf \"#MUT\taa_change\tchr\tpos\tref\talt\t${bindpred_header}\" > $out_file";

	# clean up input files for joining
	system "cat $evf | LC_ALL=C sort -k1,1 | cut -f 1,3,4,5,7,8 > ${evf}.tmp";
	system "cat $bindpred_txt | grep -v '^line_id' | LC_ALL=C sort -k1,1 > ${bindpred_txt}.tmp";
	say "$bindpred_txt";

	# join, remove seq col, add sample and mut cols, sort by consensus_percentile_rank
	my $join_cmd = 'LC_ALL=C join -t $\'\t\' -a 2';
	$join_cmd .= " ${evf}.tmp";
	$join_cmd .= " ${bindpred_txt}.tmp";
	$join_cmd .= " | cut -f 2-";
	$join_cmd .= ' | awk -F $\'\t\' \'BEGIN {OFS=FS} {print $2":"$3":"$4":"$5,$0}\'';
	$join_cmd .= " | LC_ALL=C sort -k13,13n -k14,14n";
	$join_cmd .= " >> $out_file";
	system $join_cmd;

	# delete temp files
	unlink "${evf}.tmp";
	unlink "${bindpred_txt}.tmp";

	return $out_file;
}



# end
