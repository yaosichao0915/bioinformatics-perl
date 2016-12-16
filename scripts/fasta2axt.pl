#!/usr/bin/perl

use strict;
use warnings;

use Bio::SeqIO;

my $usage = '

This script transforms fasta alignment files with 3 sequences to the axt format
required by KaKs_calculator

USAGE: perl fasta2axt.pl *.fasta

';

die $usage unless @ARGV;

while (my $fasta = shift @ARGV) {
	my @header = (); #contains the strings that are separated by _ in the fasta header. E.g. >Hel_ECU01_0460
	my @species = (); #contains the 1st string of the header. E.g. Hel
	my $gene = (); #contains the 2nd string of the header. E.g. ECU01_0460
	my @sequences = (); #contains the sequences
	my $in = Bio::SeqIO->new (-file => "$fasta", -format => 'fasta');
	while (my $seq = $in->next_seq()) {
		my $seq_names = $seq->display_id;	
		@header = split ('_', $seq_names); #split the fasta header on _
		push (@species, $header[0]);
		$gene = "$header[1]"."_"."$header[2]";
		my $sequence = $seq->seq();
		$sequence = uc($sequence);
		push (@sequences, $sequence);		
	}
	# make the pairwise files
	open OUT1, ">$gene"."_"."$species[0]"."_"."$species[1].axt";
	open OUT2, ">$gene"."_"."$species[0]"."_"."$species[2].axt";
	open OUT3, ">$gene"."_"."$species[1]"."_"."$species[2].axt";
	print OUT1 "$gene"."_"."$species[0]"."_"."$species[1]\n$sequences[0]\n$sequences[1]\n\n";
	print OUT2 "$gene"."_"."$species[0]"."_"."$species[2]\n$sequences[0]\n$sequences[2]\n\n";
	print OUT3 "$gene"."_"."$species[1]"."_"."$species[2]\n$sequences[1]\n$sequences[2]\n\n";
}
	
exit;