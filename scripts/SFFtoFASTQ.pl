#!/usr/bin/perl

use strict;
use warnings;

## Requirements: sffinfo from the gsAssembler site

my $usage = 'USAGE: SFFtoFASTQ.pl *.sff';
die $usage unless @ARGV;

## Score conversion (Sanger to Illumina, code in ASCII:Phred 33+)
my %phredQ = (
'0'=>'!','1'=>'"','2'=>'#','3'=>'$','4'=>'%','5'=>'&',
'6'=>'\'','7'=>'(','8'=>')','9'=>'*','10'=>'+','11'=>',',
'12'=>'-','13'=>'.','14'=>'/','15'=>'0','16'=>'1','17'=>'2',
'18'=>'3','19'=>'4','20'=>'5','21'=>'6','22'=>'7','23'=>'8',
'24'=>'9','25'=>':','26'=>';','27'=>'<','28'=>'=','29'=>'>',
'30'=>'?','31'=>'@','32'=>'A','33'=>'B','34'=>'C','35'=>'D',
'36'=>'E','37'=>'F','38'=>'G','39'=>'H','40'=>'I');
## End of hash

while (my $file = shift@ARGV){
        system "sffinfo $file > $file.txt";
        open IN, "<$file.txt";
        $file =~ s/.sff//;
        open OUT, ">$file.fastq";
		my $seqname = undef;   
		 my $sequence = undef;	
        while (my $line = <IN>){         
                chomp $line;
                if ($line =~ /^>(\S+)/){
                    $seqname = $1;
                    print OUT "\@$seqname\n";   ## Write down the first line on fastq format: @"name of sequence"
                }
                elsif ($line =~ /^Bases:\s+(\S+)/){
                     $sequence = $1;                    ## Extract the sequence and write it down on second line on fastq format
                    print OUT "$sequence\n";   
                    print OUT "+$seqname\n";            ## Third line is +"name of sequence"
                }                                               
                elsif ($line =~ /^Quality Scores:\s+(\d+.*$)/){ ## Extract the score of a sequence and define it in a value
                    my $allQ = $1;
                    my $quality = undef;
                    if ($allQ =~ m/(\d{1,2})(\s+)/){  ## Find and pull out the first score from the entire scores
                        $quality = $1;
                        my $key = $quality;
                        print OUT "$phredQ{$key}";      ## Score converting and output it to fastq file at the corresponding base site
                        while ($allQ =~ s/^\d{1,2}\s+//){       ## Chop the pattern of score (from left to right one by one in this loop) which has been output
                            if ($allQ =~ m/(\d{1,2})/){ ## Report the rest of scores
                                $quality = $1;
                                my $key = $quality;
                                print OUT "$phredQ{$key}";
                                }
							}
						}
						print OUT "\n";
                }
        }
       system "rm $file.sff.txt";
}

