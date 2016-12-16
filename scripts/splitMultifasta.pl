#!/usr/bin/perl

while(<>){
	BEGIN{ $/=">"; }
	if(/^\s*(\S+)/){ 
		open OUT,">$1.fsa";
		chomp;
		print OUT ">", $_
 	}
}
system "rm .fsa";