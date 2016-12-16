#!/usr/bin/perl

use strict;
use warnings;

my $usage = '

Run Muscle in a batch mode on a bunch of fasta files

USAGE: perl run_muscle.pl *.fasta

';

die $usage unless @ARGV;

while (my $fasta_file = shift @ARGV) {
	system "muscle -in $fasta_file -stable -out $fasta_file.out";
	print "aligning $fasta_file...\n";
}

exit;
