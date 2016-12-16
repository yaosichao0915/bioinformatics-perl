#!/usr/bin/perl

## Pombert Lab, IIT, 2013
## Version 1.2

use strict;
use warnings;

my $tp = 'GTGTATAAGAGACAG';		## Define transposon sequence to be removed here. Can be shorter in 5 prime (e.g. GTATAAGAGACAG), but doing so might remove good stuff. 3 prime MUST be kept. Looking only for perfect matches.
my $Ranchor = '125'; 			## Define maximum distance for the 3' end to look for TPs, Max = length of read - length of TP.
my $Lanchor = '50';			## Defined maximum distance for the 5' end to look for TPs, Max = length of read - length of TP. 

my $revtp = reverse($tp);		## Defining the reverse complement
$revtp =~ tr/ATGCatgc/TACGtacg/;	## Defining the reverse complement

my $usage = 'perl removeTP.pl *.fastq';
die $usage unless @ARGV;

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.fastq//;
	open OUT, ">$file.noTP.fastq";
	open OUT2, ">$file.noTP.stats";
	
	my $len = undef;	## set length of sequence to be conserved
	my $pos = 'good';	## set location fo sequence to keep: 5prime or 3 prime DEFAULT: good (read doesn't have TPs)
	my $seq = undef;	## sequence to be kept
	my $hed = 0;		## size of seq + tp if located at start of seq, will be set automatically, required later to set the substring offset
	my $readsnum = 0;	## initialize total number of reads
	my $readsmod = 0;	## initialize number of reads fltered
	
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^\@/){print OUT "$line\n"; $readsnum++;}	## FASTQ header
		elsif ($readsnum % 100000 == 0){print "Processing $file.fastq read number $readsnum ...\n";}
		elsif ($line =~ /^\+/){print OUT "$line\n";}	## FASTQ linker
		elsif ($line =~ /^(.*)$revtp.{0,$Ranchor}$/){	## Looking for TP in 3 prime
			$seq = $1;
			$len = length($seq);
			print OUT "$seq\n";
			$pos = '5prime';
			$readsmod++;
		}
		elsif ($line =~ /^(.{0,$Lanchor}$tp)(.*)$/){	## Looking for TP in 5 prime
			$hed = $1;
			$seq = $2;
			$len = length($seq);
			print OUT "$seq\n";
			$pos = '3prime';
			$readsmod++;
		}
		elsif ($pos eq '5prime'){			## Keeping the corresponding quality scores (if TP in 3 prime)
			my $offset = ($len); 
			my $good = substr($line, 0, $offset);
			print OUT "$good\n";
			$len = undef;
			$pos = 'good';
		}
		elsif ($pos eq '3prime'){			## Keeping the corresponding quality scores (if TP in 5 prime)
			my $offset = length($hed);
			my $good = substr($line, $offset);
			print OUT "$good\n";
			$len = undef;
			$pos = 'good';
		}
		else{print OUT "$line\n";}			## Print if no conditions are met (perfect sequences and their QS)
	}
	my $percentage = (($readsmod/$readsnum)*100);
	my $percent = sprintf("%.2f", $percentage);
	print OUT2 "Total number of reads processed from $file.fastq: $readsnum\n";
	print "Total number of reads processed from $file.fastq: $readsnum\n";
	print OUT2 "Total number of reads filtered from $file.fastq: $readsmod\n";
	print "Total number of reads filtered from $file.fastq: $readsmod\n";
	print OUT2 "Percentage of reads from the dataset $file.fastq that contained the transposon: $percent\%\n";
	print "Percentage of reads from the dataset $file.fastq that contained the transposon: $percent\%\n";
	close IN;
	close OUT;
}
exit;
