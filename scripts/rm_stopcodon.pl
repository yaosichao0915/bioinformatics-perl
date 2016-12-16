#!/usr/bin/perl

use strict;
use warnings;

use Bio::SeqIO;

my $usage = '

Cut the last 3 nucleotides in fasta files

USAGE: perl rm_stopcodon.pl fasta_list

';

my $fasta_list = shift;

open IN, "<$fasta_list" or die "can't open $fasta_list";
while (my $fasta = <IN>) {
	chomp $fasta;
	my $in = Bio::SeqIO -> new (-file => "$fasta", -format => 'fasta');
	open OUT , ">$fasta.cut";
	while (my $seq = $in -> next_seq()) {
		my $id = $seq -> id();
		my $length = $seq -> length;
		print "Gene beeing looked at: $fasta\n";
		print ">$id\n";
		print "total length = $length\n";
		my $lastbase = ($length-3);		
		print "$id has a new length of $lastbase\n";
		my $shortseq = $seq -> subseq(1,$lastbase), "\n";
		print OUT ">$id\n";
		print OUT "$shortseq\n";
	}
}

exit;