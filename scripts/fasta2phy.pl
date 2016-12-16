#!/usr/bin/perl

my $usage = 'fasta2phy.pl *.fasta';

die $usage unless @ARGV;

while (my $file =shift@ARGV){
	open FILE, "<$file";
	$file =~ s/.fasta//;
	$file =~ s/.fsa//;
	open OUT, ">$file.phy";
	
	my $taxa = -1;

	while(my $line = <FILE>){
		if($line =~ />/){    
			$taxa++;
			my $name = $line;
			$name =~ s/[\s+|\(|\)|\,|;]//g;
			$name =~ s/,//g;
			$name =~ s/>//g;
			$taxonNames[$taxa] = $name;
		}
		else{
			my $seq = $line;
			$seq =~ s/\s+//g;
			$sequences[$taxa] = $sequences[$taxa].$seq;
		}      
	}

	for($i = 0; $i <= $taxa; $i++){
		print $taxonNames[$i]." ".(length($sequences[$i]))."\n";
	}

	my $s  = $taxa + 1;
	my $bp = length($sequences[0]);
	print OUT "$s".' '."$bp\n";

	for($i = 0; $i <= $taxa; $i++){
	print OUT $taxonNames[$i]." ".$sequences[$i]."\n";
	}
}