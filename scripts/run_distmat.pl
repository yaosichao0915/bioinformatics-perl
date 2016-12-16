#!/usr/bin/perl

use strict;
use warnings;

my $usage = '

This script requires the EMBOSS package to be installed:
http://emboss.sourceforge.net/download/

distmat creates a distance matrix from a multiple sequence alignment.
This script works in a batch mode.

From distmat help:
Standard (Mandatory) qualifiers (* if not always prompted):
  [-sequence]          seqset     File containing a sequence alignment.
*  -nucmethod          menu       [0] Multiple substitution correction methods
                                  for nucleotides. (Values: 0 (Uncorrected);
                                  1 (Jukes-Cantor); 2 (Kimura); 3 (Tamura); 4
                                  (Tajima-Nei); 5 (Jin-Nei Gamma))
*  -protmethod         menu       [0] Multiple substitution correction methods
                                  for proteins. (Values: 0 (Uncorrected); 1
                                  (Jukes-Cantor); 2 (Kimura Protein))
  [-outfile]           outfile    [*.distmat] Phylip distance matrix output
                                  file

USAGE: perl run_distmat.pl *.out

';

die $usage unless @ARGV;

print "Provide a multiple substitution correction method\nOPTION for nucleotide:\n\t0 (Uncorrected)\n\t1 (Jukes-Cantor)\n\t2 (Kimura)\n\t3 (Tamura)\n\t4 (Tajima-Nei)\n\t5 (Jin-Nei Gamma))\n";
my $option = <STDIN>;
chomp $option;

while (my $aln_file = shift @ARGV) {
	#count the number of sequences in $file
	my $count = `grep -c \">\" $aln_file`; 
	chomp $count;
	if ($count > 2) {
		system "distmat $aln_file -nucmethod $option -outfile $aln_file.dist";
		print "$aln_file\n";
	}
	else {
		print "$aln_file contains only 1 sequence -> no distance matrix calculated\n";
		next;
	}
}

exit;
