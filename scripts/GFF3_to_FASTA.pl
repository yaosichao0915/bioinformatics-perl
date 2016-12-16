#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perl GFF3_to_FASTA.pl *.gff';

die $usage unless @ARGV;

## Convert the GFF3 output from the Maker gauntlet to the fasta (.fsa) format required by TBL2ASN ##

while (my $gff = shift @ARGV) {

	open IN, "<$gff" or die "cannot open $gff";
	$gff =~ s/\.gff$//;
	open OUT, ">$gff.fsa";

	while (my $line = <IN>) {

		chomp $line;

		if ($line =~ /^##.*$/) {
			next;
		}
		
		elsif ($line =~ /^contig-\d+.*$/) {
			next;
		}
		
		elsif ($line =~ /^>contig-\d+.*$/) {
			print OUT "$line.fsa\n";
		}
			
		else {
			print OUT "$line\n";
		}
	}
	close IN;
}

exit;
