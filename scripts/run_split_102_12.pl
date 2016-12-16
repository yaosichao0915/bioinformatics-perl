#!/usr/bin/perl

use strict;
use warnings;

my $usage = '

perl run_split.pl *.axt

';

die $usage unless @ARGV;

while (my $axt_file = shift @ARGV) {
	### java split file.axt [window size] [step size]
	system "java split $axt_file 102 12";
}

exit;
