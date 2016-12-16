#!usr/local/bin/perl

### This script takes 2 inputs: the 454ReadStatus file and the reads file in fasta format.
### It creates a fasta file named singletons.fasta containing the reads NOT assembled into contigs
### during the Newbler assembling process and that are larger than a threshold given on the command line

### USAGE: perl GetSingleton.pl 454ReadStatus.txt reads.fasta

use strict;
use warnings;
use Bio::SeqIO;

my $ReadStatus = shift or die "no 454ReadStatus.txt file provided\n";
my $FastaReads = shift or die "no reads fasta file provided\n";
my %ReadId_status = ();

print "Minimum length of a read to be retrieved:\n";
my $length = <STDIN>;

my $in  = Bio::SeqIO->new(-file => "$FastaReads" , -format => 'Fasta');
my $out = Bio::SeqIO->new(-format => 'Fasta', -file => ">singletons_$length.fasta");

open IN1, "<$ReadStatus" or die "cannot open $ReadStatus";
while (my $line = <IN1>) {
	if ($line =~ /^Accno\t/) {
		next;
	}
	elsif ($line =~ /^(\w+)\t\Singleton$/) {
		$ReadId_status{$1}=1;
	}
}

close IN1;

foreach my $key (keys %ReadId_status) {
	print "$key\n";
}

my $count = keys %ReadId_status;
print "\nFound $count Singletons\n";

my $counter = 0;
##Go through each sequence in the fasta file
while (my $seq = $in->next_seq()) {
	my $id = $seq->id();
	if (($ReadId_status{$id}) and ($seq->length() >= $length)) {
		$out->write_seq($seq); 
		$counter = $counter + 1;
	}
}

print "Found $counter reads >= $length\n\n";

exit;