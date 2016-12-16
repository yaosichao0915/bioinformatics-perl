#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'Run Mafft in a batch mode on a bunch of fasta files USAGE: perl run_mafft.pl *.fasta';

die $usage unless @ARGV;

while (my $fasta = shift@ARGV) {
	print "aligning $fasta...\n";
	system "mafft --thread 32 $fasta > $fasta.aln";
}

exit;
