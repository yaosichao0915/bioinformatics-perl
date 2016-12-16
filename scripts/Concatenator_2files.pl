#!/usr/bin/perl

### Loop through 2 fasta files and generates files containing sequences found in all parent files
### !!! It is hard coded to match a given fasta header, so won't work in all instances !!!
### This version works for header like this: >Ecu01_0025-0030 8833:11157 forward
### The tag that is searched and must be present in all parent files is 0025-0030

### USAGE: perl Concatenator.pl fasta_file1 fasta_file2

use strict;
use warnings;

use Bio::SeqIO;
$ENV{BIOPERL_INDEX_TYPE} = "SDBM_File";
use Bio::Index::Fasta;
use List::Compare;
use Data::Dumper;

my $input1 = shift;
my $input2 = shift;

#make index of the fasta files for fast retrieving, if already exist skip it
my $idx1 = Bio::Index::Fasta->new(-filename => $input1 . ".idx", -write_flag => 1);
$idx1->id_parser(\&get_id);
$idx1->make_index($input1);
my $idx2 = Bio::Index::Fasta->new(-filename => $input2 . ".idx", -write_flag => 1);
$idx2->id_parser(\&get_id);
$idx2->make_index($input2);

my $seq_in1  = Bio::SeqIO->new( -format => 'fasta', -file => $input1);
my $seq_in2  = Bio::SeqIO->new( -format => 'fasta', -file => $input2);

# create arrays containing the tags that need to be identical throughout the files
my (@tags1, @tags2);
while(my $seq_obj1 = $seq_in1->next_seq() ) {
	my $id1 = $seq_obj1->display_id;
	if ($id1 =~ /\S+_(\S+)/) {
		my $tag1 = $1;
		push(@tags1,$tag1);
	}
}
while(my $seq_obj2 = $seq_in2->next_seq() ) {
	my $id2 = $seq_obj2->display_id;
	if ($id2 =~ /\S+_(\S+)/) {
		my $tag2 = $1;
		push(@tags2,$tag2);
	}
}

# use List::Compare to find common tags across the fasta files
my $lc = List::Compare->new(\@tags1, \@tags2);
my @common = $lc->get_intersection;

# fetch the sequences and create files according to the list of common tags
foreach my $common (@common) {
	my $out1 = Bio::SeqIO->new(-file => ">>$common.fasta", -format => 'Fasta');
	my $out2 = Bio::SeqIO->new(-file => ">>$common.fasta", -format => 'Fasta');
	my $seq_obj1 = $idx1->fetch($common);
	my $seq_obj2 = $idx2->fetch($common);
	$out1->write_seq($seq_obj1);
	$out2->write_seq($seq_obj2);
}

### subroutine section ###

sub get_id {       
   	my $header = shift;       
  	$header =~ /^>\S+_(\S+)\s.*$/;       
    $1;
}
	