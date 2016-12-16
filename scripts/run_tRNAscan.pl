#!/usr/bin/perl

use strict;
use warnings;

while (my $file = shift@ARGV){
	system "tRNAscan-SE $file > $file.tRNAs";
}