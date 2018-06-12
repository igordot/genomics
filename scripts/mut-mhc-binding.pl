#!/usr/bin/env perl

use strict;
use warnings;
use v5.10;
use File::Basename;
use Bio::SeqIO;

my $HELP = <<HELP;

  Predict MHC-I and MHC-II binding for ANNOVAR-annotated point mutations using IEDB MHC binding prediction tool.
  Protein coding changes are annotated using ANNOVAR.
  Binding predictions are calculated by IEDB MHC class I/II binding prediction tool using consensus method.
  MHC class is selected automatically based on the binding tool.

  Usage:
    perl mut-mhc-binding.pl base_name avinput_file annovar_scripts_dir annovar_refgene_txt iedb_mhc_dir

  Arguments:
    base_name             output prefix
    avinput_file          ANNOVAR avinput file - 5 required cols (chr, start, end, ref, alt) and then optional cols
    annovar_scripts_dir   ANNOVAR scripts directory
    annovar_refgene_txt   ANNOVAR refGene.txt reference file
    iedb_mhc_dir          IEDB MHC I or II binding tool directory (tested with version 2.17)

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
	my $iedb_mhc_dir = $ARGV[4];

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
	unless ( -d $iedb_mhc_dir ) {
		die "\n\n ERROR: $iedb_mhc_dir does not exist \n\n";
	}

	# check that the ANNOVAR scripts exist
	my $annovar_annotate_pl = "${annovar_scripts_dir}/annotate_variation.pl";
	unless ( -e $annovar_annotate_pl ) {
		die "\n\n ERROR: $annovar_annotate_pl does not exist \n\n";
	}
	my $annovar_coding_change_pl = "${annovar_scripts_dir}/coding_change.pl";
	unless ( -e $annovar_coding_change_pl ) {
		die "\n\n ERROR: $annovar_coding_change_pl does not exist \n\n";
	}

	# peptide length for binding predictions (and preparing binding prediction input sequences)
	my $peptide_length;

	# check that the IEDB scripts exist and determine class based on the available script
	my $mhc_class;
	my $predict_binding_py;
	my $predict_binding_i_py = "${iedb_mhc_dir}/src/predict_binding.py";
	my $predict_binding_ii_py = "${iedb_mhc_dir}/mhc_II_binding.py";
	if ( -e $predict_binding_i_py ) {
		$mhc_class = "I";
		$predict_binding_py = $predict_binding_i_py;
		# 8-14 on the IEDB and NetMHC web server
		$peptide_length = 9;
	}
	elsif ( -e $predict_binding_ii_py ) {
		$mhc_class = "II";
		$predict_binding_py = $predict_binding_ii_py;
		# generally between 15 and 24, 15 is default on the NetMHC web server
		$peptide_length = 15;
	}
	else {
		die "\n\n ERROR: binding.py does not exist \n\n";
	}

	# update file prefix since the rest will depend on MHC class
	$base_name = "${base_name}.MHC-${mhc_class}.pad${peptide_length}";

	say "annotating input mutations (avinput file)";
	my $evf = annotate_avinput($base_name, $avinput, $annovar_annotate_pl, $annovar_refgene_txt);

	say "getting coding change";
	my $cod_change_fa = get_coding_change($base_name, $evf, $annovar_coding_change_pl, $annovar_refgene_txt);

	say "formatting coding change";
	my $mutpad_fa = format_coding_change_fa($base_name, $cod_change_fa, $peptide_length);

	say "predicting binding";
	my $bindpred_txt;
	if ($mhc_class eq "I") {
		say "MHC class: I";
		$bindpred_txt = predict_mhc_i_binding($base_name, $mutpad_fa, $peptide_length, $predict_binding_py);
	}
	if ($mhc_class eq "II") {
		say "MHC class: II";
		$bindpred_txt = predict_mhc_ii_binding($base_name, $mutpad_fa, $peptide_length, $predict_binding_py);
	}

	say "annotating binding predictions";
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
	if ( -z $evf ) {
		die "\n\n ERROR: $evf IS EMPTY \n\n";
	}

	# clean up
	unlink "${base_name}.log";
	unlink "${base_name}.variant_function";

	say "created $evf";

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
	if ( -z $out_file ) {
		die "\n\n ERROR: $out_file IS EMPTY \n\n";
	}

	say "created $out_file";

	return $out_file;
}

# format coding change FASTA for IEDB binding predictions
sub format_coding_change_fa {
	my $base_name = $_[0];
	my $coding_change_fa = $_[1];
	my $peptide_length = $_[2];

	# sequence padding (peptide length)
	my $padding = $peptide_length - 1;

	my $out_file = "${base_name}.fa";

	# delete output if already exists
	if ( -e $out_file ) {
		unlink $out_file;
	}

	my $seqio_in = Bio::SeqIO->new(-file => $coding_change_fa, -format => 'fasta');
	my $seqio_out = Bio::SeqIO->new(-file => ">$out_file", -format => 'fasta');
	my $seq_num = 0;
	# tracking WT sequence info (listed before mutant)
	my ($seq_id_wt, $tx_name_wt, $sequence_wt, $seq_padded_wt) = ("ERR", "ERR", "ERR", "ERR");
	while ( my $seq_in = $seqio_in->next_seq() ) {
		# extract relevant parts
		my $seq_id = $seq_in->id;
		my $seq_desc = $seq_in->desc;
		$seq_desc =~ s/\s+$//;
		my $sequence = $seq_in->seq;

		# extract mutation details and sequence around the WT sequence (listed before mutant)
		if ($seq_desc =~ m/(.+?)\sWILDTYPE/) {
			$seq_id_wt = $seq_id;
			$sequence_wt = $sequence;
			$tx_name_wt = $1;
		}

		# process only altered sequences and single amino acid substitutions
		if ($seq_desc !~ m/silent/ && $seq_desc =~ m/\sp\.\w\d+\w\s/) {
			my ($tx_name, $aa_num, $aa_from, $aa_to, $seq_padded) = ("ERR", "ERR", "ERR", "ERR", "ERR");
			# extract mutation details and sequence around the altered amino acid
			if ($seq_desc =~ m/(.+?)\s.+?\(position\s(\d+)\schanged\sfrom\s(\w)\sto\s(\w)\)/) {
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

				# repeat for WT (was processed as previous sequence)
				my $seq_padded_wt = substr($sequence_wt, $substr_start, $substr_length);
				$seq_padded_wt =~ s/\*$//;

				# create output FASTA record and write it (binding predictions will just assign consecutive numbers)
				$seq_num++;
				my $new_id = "${seq_num}|${seq_id}|${tx_name}|${aa_from}${aa_num}${aa_to}|${seq_padded}";
				my $seq_out = Bio::Seq->new(-seq => $seq_padded, -id => $new_id);
				$seqio_out->write_seq($seq_out);

				# create WT output FASTA record and write it (binding predictions will just assign consecutive numbers)
				$seq_num++;
				my $new_id_wt = "${seq_num}|${seq_id_wt}|${tx_name_wt}|WILDTYPE|${seq_padded_wt}";
				my $seq_out_wt = Bio::Seq->new(-seq => $seq_padded_wt, -id => $new_id_wt);
				$seqio_out->write_seq($seq_out_wt);
			}
		}
	}

	# confirm that padded mutations FASTA file generated
	unless ( -e $out_file ) {
		die "\n\n ERROR: $out_file DOES NOT EXIST \n\n";
	}
	if ( -z $out_file ) {
		die "\n\n ERROR: $out_file IS EMPTY \n\n";
	}

	say "created $out_file";

	return $out_file;
}

# run IEDB MHC-I Binding Predictions
sub predict_mhc_i_binding {
	my $base_name = $_[0];
	my $mutpad_fa = $_[1];
	my $peptide_length = $_[2];
	my $predict_binding_py = $_[3];

	my $raw_iedb_out_file = "${base_name}.iedb.txt";

	# delete output if already exists
	if ( -e $raw_iedb_out_file ) {
		unlink $raw_iedb_out_file;
	}

	# a reference panel of 27 alleles (human HLA reference set with maximal population coverage)
	# http://help.iedb.org/hc/en-us/articles/114094151851
	# my @alleles = ('HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*02:03', 'HLA-A*02:06', 'HLA-A*03:01', 'HLA-A*11:01',
	#	'HLA-A*23:01', 'HLA-A*24:02', 'HLA-A*26:01', 'HLA-A*30:01', 'HLA-A*30:02', 'HLA-A*31:01','HLA-A*32:01',
	#	'HLA-A*33:01', 'HLA-A*68:01', 'HLA-A*68:02', 'HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*15:01', 'HLA-B*35:01',
	#	'HLA-B*40:01', 'HLA-B*44:02', 'HLA-B*44:03', 'HLA-B*51:01', 'HLA-B*53:01', 'HLA-B*57:01', 'HLA-B*58:01');

	# mouse alleles
	my @alleles = ('H-2-Db', 'H-2-Dd', 'H-2-Kb', 'H-2-Kd', 'H-2-Kk', 'H-2-Ld');

	# make predictions for each allele
	foreach (@alleles) {
		say "running predictions for allele $_";

		# predict_binding.py command
		# ./src/predict_binding [method] [mhc] [peptide_length] [input_file]
		my $predict_binding_cmd = "$predict_binding_py consensus $_ $peptide_length $mutpad_fa >> $raw_iedb_out_file";
		system $predict_binding_cmd;

		# confirm that binding predictions file generated
		unless ( -e $raw_iedb_out_file ) {
			die "\n\n ERROR: $raw_iedb_out_file DOES NOT EXIST \n\n";
		}
		if ( -z $raw_iedb_out_file ) {
			die "\n\n ERROR: $raw_iedb_out_file IS EMPTY \n\n";
		}
	}

	my $merged_output = merge_binding_predictions_input_output($base_name, $mutpad_fa, $raw_iedb_out_file);

	# clean up
	sleep(1);
	# unlink $mutpad_fa;

	say "created $merged_output";

	return $merged_output;

}

# run IEDB MHC-II Binding Predictions
sub predict_mhc_ii_binding {
	my $base_name = $_[0];
	my $mutpad_fa = $_[1];
	my $peptide_length = $_[2];
	my $predict_binding_py = $_[3];

	my $raw_iedb_out_file = "${base_name}.iedb.txt";

	# delete output if already exists
	if ( -e $raw_iedb_out_file ) {
		unlink $raw_iedb_out_file;
	}

	# mouse alleles
	my @alleles = ('H2-IAb', 'H2-IAd');

	# make predictions for each allele
	foreach (@alleles) {
		say "running predictions for allele $_";

		# mhc_II_binding.py command
		# python mhc_II_binding.py prediction_method_name allele_name input_sequence_file_name
		# Example: python mhc_II_binding.py consensus3 HLA-DRB1*03:01 test.fasta
		my $predict_binding_cmd = "$predict_binding_py consensus3 $_ $mutpad_fa >> $raw_iedb_out_file";
		system $predict_binding_cmd;
		say "$predict_binding_cmd";

		# confirm that binding predictions file generated
		unless ( -e $raw_iedb_out_file ) {
			die "\n\n ERROR: $raw_iedb_out_file DOES NOT EXIST \n\n";
		}
		if ( -z $raw_iedb_out_file ) {
			die "\n\n ERROR: $raw_iedb_out_file IS EMPTY \n\n";
		}
	}

	my $merged_output = merge_binding_predictions_input_output($base_name, $mutpad_fa, $raw_iedb_out_file);

	# clean up
	sleep(1);
	# unlink $mutpad_fa;

	say "created $merged_output";

	return $merged_output;

}

# combine binding predictions input FASTA and output table
sub merge_binding_predictions_input_output {
	my $base_name = $_[0];
	my $mutpad_fa = $_[1];
	my $iedb_out_txt = $_[2];

	my $out_file = "${base_name}.binding.txt";

	# delete output if already exists
	if ( -e $out_file ) {
		unlink $out_file;
	}

	# $iedb_out_txt columns: allele, seq_num, start, end, peptide, ...
	# ${mutpad_fa}.txt columns: seq_num, line_id, transcript_id, mutation, peptide

	# header for the mutated and WT parts of the combined file
	my $bindpred_header = `cat $iedb_out_txt | head -1 | cut -f 1,3,5-11`;
	system "printf \"line_id\ttranscript_id\tp_change\taa_padded\t${bindpred_header}\" > ${out_file}.mut.txt";
	$bindpred_header =~ s/\t/_wt\t/g;
	system "printf \"p_change_wt\taa_padded_wt\t${bindpred_header}\" > ${out_file}.wt.txt";

	# clean up input files for joining
	system "cat $mutpad_fa | grep '^>' | cut -c 2- | tr '|' '\t' | LC_ALL=C sort -k1,1 > ${mutpad_fa}.txt";
	system "cat $iedb_out_txt | grep -v '^allele' | LC_ALL=C sort -k2,2 -k1,1 -k3,3 | cut -f 1,2,3,5-11 > ${iedb_out_txt}.tmp";

	system "cat ${mutpad_fa}.txt | grep -v 'WILDTYPE' > ${mutpad_fa}.mut.txt";
	system "cat ${mutpad_fa}.txt | grep 'WILDTYPE' > ${mutpad_fa}.wt.txt";

	# join by seq_num, remove seq_num col, sort by line_id and start
	my $join_mut_cmd = 'LC_ALL=C join -t $\'\t\' -a 1 -1 1 -2 2';
	$join_mut_cmd .= " ${mutpad_fa}.mut.txt";
	$join_mut_cmd .= " ${iedb_out_txt}.tmp";
	$join_mut_cmd .= " | LC_ALL=C sort -k2,2 -k7,7";
	$join_mut_cmd .= " | cut -f 2-";
	$join_mut_cmd .= " >> ${out_file}.mut.txt";
	system $join_mut_cmd;

	# join by seq_num
	my $join_wt_cmd = 'LC_ALL=C join -t $\'\t\' -a 1 -1 1 -2 2';
	$join_wt_cmd .= " ${mutpad_fa}.wt.txt";
	$join_wt_cmd .= " ${iedb_out_txt}.tmp";
	$join_wt_cmd .= " | LC_ALL=C sort -k2,2 -k7,7";
	$join_wt_cmd .= " | cut -f 4-";
	$join_wt_cmd .= " >> ${out_file}.wt.txt";
	system $join_wt_cmd;

	sleep(1);

	system "paste ${out_file}.mut.txt ${out_file}.wt.txt >> ${out_file}";

	# confirm that binding predictions file generated
	unless ( -e $out_file ) {
		die "\n\n ERROR: $out_file DOES NOT EXIST \n\n";
	}

	# clean up
	sleep(1);
	# unlink $mutpad_fa;
	unlink "${iedb_out_txt}.tmp";
	unlink "${mutpad_fa}.mut.txt";
	unlink "${mutpad_fa}.wt.txt";
	unlink "${out_file}.mut.txt";
	unlink "${out_file}.wt.txt";

	return $out_file;

}

# combine binding predictions table with variant annotations
sub annotate_binding_predictions {
	my $base_name = $_[0];
	my $bindpred_txt = $_[1];
	my $evf = $_[2];

	my $out_file = "${base_name}.binding.annot.txt";

	# header for the combined file
	my $bindpred_header = `cat $bindpred_txt | head -1 | cut -f 2-`;
	system "printf \"#MUT\taa_change\tchr\tpos\tref\talt\t${bindpred_header}\" > $out_file";

	# clean up input files for joining
	system "cat $evf | LC_ALL=C sort -k1,1 | cut -f 1,3,4,5,7,8 > ${evf}.tmp";
	system "cat $bindpred_txt | grep -v '^line_id' | LC_ALL=C sort -k1,1 > ${bindpred_txt}.tmp";

	# join, remove seq col, add sample and mut cols, sort by consensus_percentile_rank
	my $join_cmd = 'LC_ALL=C join -t $\'\t\' -a 2';
	$join_cmd .= " ${evf}.tmp";
	$join_cmd .= " ${bindpred_txt}.tmp";
	$join_cmd .= " | cut -f 2-";
	$join_cmd .= ' | awk -F $\'\t\' \'BEGIN {OFS=FS} {print $2":"$3":"$4":"$5,$0}\'';
	$join_cmd .= " | LC_ALL=C sort -k13,13n -k14,14n";
	$join_cmd .= " >> $out_file";
	system $join_cmd;

	# clean up
	sleep(1);
	unlink "${evf}.tmp";
	unlink "${bindpred_txt}.tmp";

	return $out_file;
}



# end
