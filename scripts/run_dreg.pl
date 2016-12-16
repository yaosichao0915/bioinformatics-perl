#!/usr/bin/perl

use strict;
use warnings;

my $usage = '

This script requires the EMBOSS package to be installed:
http://emboss.sourceforge.net/download/

Regular expression search of nucleotide sequence(s).
This script works in a batch mode.

>From dreg help:

Standard (Mandatory) qualifiers:
	[-sequence]          seqall     Nucleotide sequence(s) filename and optional
                                  format, or reference (input USA)
	[-pattern]           regexp     Any regular expression pattern is accepted)
	[-outfile]           report     [*.dreg] Output report file name
  
Additional (Optional) qualifiers (among many):
	-sreverse1          boolean    Reverse (if DNA)
  
USAGE: perl run_dreg.pl *.fasta

';

die $usage unless @ARGV;

while (my $sequence_file = shift @ARGV) {
	system "dreg -sequence $sequence_file -pattern 'gta[ag]gt[acgt]\{5,30\}tt[acgt]\{0,3\}AG' -outfile $sequence_file".".dreg";
	system "dreg -sequence $sequence_file -sreverse1 -pattern 'gta[ag]gt[acgt]\{5,30\}tt[acgt]\{0,3\}AG' -outfile $sequence_file"."_reverse.dreg";
}

exit;
