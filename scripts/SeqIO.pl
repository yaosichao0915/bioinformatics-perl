#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

## Based on BioPerl's SeqIO module. Requires local installation of BioPerl
## See www.bioperl.org/wiki/HOWTO:SeqIO for supported formats, more information and credits
## Common formats supported are: fasta, fastq, abi, ace, embl, genbank, interpro, kegg, excel, tab, table

my $usage = 'USAGE: SeqIO.pl input_format output_format *.files';

my $format = shift@ARGV or die "usage = $usage\n"; ## input file format
my $outfmt = shift@ARGV or die "usage = $usage\n"; ## output file format

while (my $file = shift@ARGV){

	my $in = Bio::SeqIO->new(-file => "$file", -format => "$format");
	$file =~ s/\.$format//;
	my $out = Bio::SeqIO->new(-file => ">$file.$outfmt", -format => "$outfmt");
	
	while (my $seq = $in->next_seq()){
		$out->write_seq($seq);
	}
}
