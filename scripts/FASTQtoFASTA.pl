#!/usr/bin/perl
## Fan Zhao, Pombert Lab, IIT 2014

use strict;
use warnings;

my $usage = 'perl FASTQ_to_FASTA *.fastq';
die $usage unless @ARGV;

### Quality scores hashes 33/64 conversion
my %sanger33 = ('!'=>'0','"'=>'1','#'=>'2','$'=>'3','%'=>'4','&'=>'5','\''=>'6','('=>'7',')'=>'8','*'=>'9','+'=>'10',','=>'11','-'=>'12','.'=>'13','/'=>'14','0'=>'15','1'=>'16','2'=>'17','3'=>'18','4'=>'19','5'=>'20','6'=>'21','7'=>'22','8'=>'23','9'=>'24',':'=>'25',';'=>'26','<'=>'27','='=>'28','>'=>'29','?'=>'30','@'=>'31','A'=>'32','B'=>'33','C'=>'34','D'=>'35','E'=>'36','F'=>'37','G'=>'38','H'=>'39','I'=>'40',);
my %illumina64 = ('@'=>'0','A'=>'1','B'=>'2','C'=>'3','D'=>'4','E'=>'5','F'=>'6','G'=>'7','H'=>'8','I'=>'9','J'=>'10','K'=>'11','L'=>'12','M'=>'13','N'=>'14','O'=>'15','P'=>'16','Q'=>'17','R'=>'18','S'=>'19','T'=>'20','U'=>'21','V'=>'22','W'=>'23','X'=>'24','Y'=>'25','Z'=>'26','['=>'27','\\'=>'28',']'=>'29','^'=>'30','_'=>'31','`'=>'32','a'=>'33','b'=>'34','c'=>'35','d'=>'36','e'=>'37','f'=>'38','g'=>'39','h'=>'40',);
### End of hashes

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.fastq//;
	open OUT1, ">$file.fasta";
	open OUT2, ">$file.qual";

	my $readline = 1;
	my $score_type = 0;	## Undef;
	my @Quality = ();

	while (my $line = <IN>){	## Output line by line to fasta format in one set
		chomp $line;
		if ($readline == 1){	## Output name of sequence
		    $line =~ s/@/>/;
		    print OUT1 "$line\n";
		    print OUT2 "$line\n";
		    $readline++;
		}
		elsif($readline == 2){          ## Output sequence
		      print OUT1 "$line\n";
		      $readline++;
		}
		elsif($readline == 3){
		      $readline++;
		      next;
		}
		elsif(($readline == 4)&&($score_type == 0)){
		      	 @Quality = split('', $line);			## Crack the quality line to individual TIMTOWDI "unpack ("(A1)*", $line)"
				while (my $search = shift(@Quality)){	## One by one character matching search till found, to determine type
					if ((exists $sanger33{$search})&&(exists $illumina64{$search})){
					next;	## We don't know
					}
					elsif (exists $sanger33{$search}){   ## Detect first set of data meanwhile output it
					    $score_type = 33;
 					    @Quality = split('', $line);   
							foreach my $Quality (@Quality){
                                my $key_33 = $Quality;
                                print OUT2 "$sanger33{$key_33}\t";
                            }
                        print OUT2 "\n";
		    			$readline = 1;
					last;                                
					}                                    ## End of 33 detect
					elsif (exists $illumina64{$search}){ ## Detect first set of data meanwhile have to output it, otherwise readline messed up
					    $score_type = 64;
						@Quality = split('', $line); 
							foreach my $Quality (@Quality){
								my $key_64 = $Quality;
								print OUT2 "$illumina64{$key_64}\t";
							}
						print OUT2 "\n";
					    $readline = 1;
					last;
					}									## End of 64 detect
 					else {                              ## Overlap score
 					   my $usage = "Can't determine the type of quality";
 					   die $usage; 
 					}				
				}
		}
		elsif(($readline == 4)&&($score_type == 33)){   ## Sanger33 to PhredQ converting
		      	 @Quality = split('', $line); 			## Crack the quality line to individual TIMTOWDI "unpack ("(A1)*", $line)"
				foreach my $Quality (@Quality){
					my $key_33 = $Quality;         
		  			print OUT2 "$sanger33{$key_33}\t";
				}
				print OUT2 "\n";
				$readline = 1;              ## Start over next one set of data				
		}									## End of 33 transfer
		elsif(($readline == 4)&&($score_type == 64)){   ## Illumina to PhredQ converting
		       @Quality = split('', $line); 			## Crack the quality line to individual TIMTOWDI "unpack ("(A1)*", $line)"
				foreach my $Quality (@Quality){
					my $key_64 = $Quality;
					print OUT2 "$illumina64{$key_64}\t";
				}
				print OUT2 "\n";
				$readline = 1;
		}	     							## End of 64 transfer
	}
print "the type of quality score is $score_type\n";
}

