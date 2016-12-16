#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perk run_rnahammer.pl kingdom files.fsa'; ## kingdom = arc, bac or euk

my $kingdom = shift@ARGV;

while (my $file = shift@ARGV){
	system "rnammer -S $kingdom -m tsu,ssu,lsu -gff $file.gff -h $file.hmm -f $file.rRNAs < $file";
}
