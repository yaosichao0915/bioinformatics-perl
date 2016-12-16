#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'USAGE: perl SFFtoFASTA.pl *.sff';
die $usage unless @ARGV;

while (my $file = shift@ARGV){

	system "sffinfo $file > $file.txt"; ## We generate a text file from the binary SFF file.

	open IN, "<$file.txt";	## We open the text file
	$file =~ s/.sff//;
	open OUT, ">$file.fasta";
	open OUT2, ">$file.qual";
 
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>(\S+)/){ ## Looking for the sequence names
			my $seqname = $1;
			print OUT ">$seqname\n";
			print OUT2 ">$seqname\n";
		}
		elsif ($line =~ /^Bases:	(\S+)/){  ## Printing sequences to $file.fasta
			my $sequence = $1;
			print OUT "$sequence\n";
		}
		elsif ($line =~ /^Quality Scores:\s+(\d+.*$)/){  ## Printing quality values to $file.qual
			my $qvs = $1;
			print OUT2 "$qvs\n";
		}
	}
	system "rm $file.sff.txt"; ## We remove the text file we created for the SFF to Fasta conversion
}
