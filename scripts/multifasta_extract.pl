#!/usr/bin/perl

use Bio::DB::Fasta;

## usage: perl thisScript.pl fastaFile.fa queryFile.txt ##

my $fastaFile = shift;
my $queryFile = shift;

my $db = Bio::DB::Fasta->new( $fastaFile );
open (IN, $queryFile);
while (<IN>){
    chomp; 
    $seq = $_;
    my $sequence = $db->seq($seq);
    if  (!defined( $sequence )) {
            die "Sequence $seq not found. \n" 
    }   
    print ">$seq\n", "$sequence\n";
}

## The first time you use it, it will index the fasta file, but then it will be superfast and will let you fetch all sequences (one ID per line in file ) ##
## The good thing is that after indexing, it KNOWS where every sequence is in the file and it will get it without going through the full file. ##
