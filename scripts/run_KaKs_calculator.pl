#!/usr/bin/perl

use strict;
use warnings;

my $usage = '

perl run_KaKs_calculator.pl *.axt

';

die $usage unless @ARGV;

while (my $axt_file = shift @ARGV) {
	system "KaKs_Calculator -i $axt_file -o $axt_file.kaks -m GMYN";
}

exit;
