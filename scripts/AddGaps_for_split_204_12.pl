#!/usr/bin/perl

use strict;
use warnings;

use Bio::SeqIO;

my $usage = '

This script adds gap characters to sequences in a fasta file so that the program SPLIT
does not truncate the last window

USAGE: perl AddGaps_for_split *.fasta

';

die $usage unless @ARGV;

while (my $fasta = shift @ARGV) {
	my $in = Bio::SeqIO->new (-file => "$fasta", -format => 'fasta');
	
	### VERIFY THAT THE WINDOW SIZE AND STEP PARAMETERS ARE CORRECT ###
	my $windows_size = '204';
	my $step = '12';
	###################################################################
	
	open OUT, ">>AddGaps.out"; #contains info on the successive try outs  until a window/step prefectly fits
	open OUT1, ">$fasta"."_gap"."."."$windows_size";
	
	print OUT "\nwindows size: $windows_size\n";
	print OUT "step: $step\n##########\n";	
	my $gap_added = ();
	while (my $seq = $in->next_seq()) {
		my @gaps = (); #contains the gaps that need to be added
		my $seq_names = $seq->display_id;
		print OUT "\n***$seq_names***\n";
		my $initial_gene_length = $seq->length;
		my $gene_length = $initial_gene_length;
		my $windows_nb = ($initial_gene_length-$windows_size)/$step; #this is the number of windows, given the gene length, the window size and the step value
		print OUT "initial gene length: $initial_gene_length\n";
		print OUT "initial window numbers: $windows_nb\n\n";
		until ($windows_nb =~ /^[0-9]{1,}$/) { #check for integer, if the number of window is not an integer, don't stop
			$gene_length=$gene_length+3; #add 3 gaps
			$gap_added=$gene_length-$initial_gene_length;
			print OUT "gap added: $gap_added\n";
			$windows_nb = ($gene_length-$windows_size)/$step;
			print OUT "new windows number: $windows_nb\n";
			print OUT "gene size with gap added: $gene_length\n";
		}
		for (my $i = 1; $i <= $gap_added; $i++) {
			push (@gaps, '-'); #make an array with the right number of gaps
		}
		### print out the sequence+gaps
		my $sequence = $seq->seq();
		$sequence = uc($sequence);
		if (($windows_size-$gap_added) >= 21) { #add gaps only if at least 21 nucleotides are present in the last window
			print OUT1 ">$seq_names\n$sequence";
			print OUT1 @gaps,"\n";
			print OUT "-----\n";
		}
		else {
			print OUT1 ">$seq_names\n$sequence\n";
			print OUT "no gap added for $seq_names, <20 nucleotides would be present in the last windows\n-----\n";
		}
	}
}

exit;
