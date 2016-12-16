#!/usr/bin/perl

use strict;
use warnings;

my $usage = '

This script parses the output of distmat

USAGE: perl DistParser.pl *.dist

';

die $usage unless @ARGV;

open (OUT, ">distances.txt");

print OUT "gene name\tseq1 vs seq2\tseq1 vs seq3\tseq2 vs seq3\n";

while (my $dist_file = shift @ARGV) {	
	print OUT "$dist_file\t";
	open (IN, "<$dist_file") or die ("can't open $dist_file");
	while (my $line = <IN>) {
		if ($line =~ /\s+0.00\s+(\S+)\s+(\S+)\s+\S+\s\d{1,}$/) {
			print OUT "$1\t$2";
		}
		elsif ($line =~ /\s+0.00\s+(\S+)\s+\S+\s\d{1,}$/) {
			print OUT "\t$1\n";
		}
	}
}

exit;
