#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perl run_macse.pl *.fasta';

die $usage unless @ARGV;

while (my $fasta_file = shift @ARGV) {
	$fasta_file =~ s/.fasta//;
	system "java -jar /opt/MACSE/macse_v0.9b1.jar -i $fasta_file.fasta -o $fasta_file.macse";
	print "aligning $fasta_file...\n";
}

exit;
