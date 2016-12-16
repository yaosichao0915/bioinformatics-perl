#!/usr/bin/perl -w

use strict;
use warnings;

my $usage = '

USAGE: perl rename_in_fasta.pl *.fasta

';

while (my $fasta = shift @ARGV) {
	open IN, "<$fasta" or die "cannot open $fasta";
	open OUT, ">$fasta.txt";
	while (my $line = <IN>) {
		chomp $line;
		if ($line =~ /^>Hel\s(.*)$/) {
			print OUT ">Hel_$1\n";
		}
		elsif ($line =~ /^>Ein\s(.*)$/) {
			print OUT ">Ein_$1\n";
		}
		elsif ($line =~ /^>Ecu\s(.*)$/) {
			print OUT ">Ecu_$1\n";
		}
		elsif ($line =~ /^>ECU(.*)$/) {
			print OUT ">Ecu_ECU$1\n";
		}
		elsif ($line =~ /^>(\w{3}\_\S+\_\d{4})(\_\d{1,})$/) {
			print OUT ">$1\n";
		}
		else {
			print OUT "$line\n";
		}
	}
	close IN;
}

exit;
