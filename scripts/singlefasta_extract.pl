#!/usr/bin/perl

use warnings;
use strict;

## usage: fasta_extract.pl contig < multifasta.fsa > output.fsa ##

my $lookup = shift @ARGV;  # ID to extract

local $/ = "\n>";  # read by FASTA record

while (my $seq = <>) {
    chomp $seq;
    my ($id) = $seq =~ /^>*(\S+)/;  # parse ID as first word in FASTA header
    if ($id eq $lookup) {
        $seq =~ s/^>*.+\n//;  # remove FASTA header
        $seq =~ s/\n//g;  # remove endlines
        print "$seq\n";
        last;
    }
}

