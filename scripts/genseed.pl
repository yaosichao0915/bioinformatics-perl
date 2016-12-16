#!/usr/bin/perl
my $version = "1.0.22 - 16-SEP-2008";

use strict;
use Getopt::Long;
use GD;

###Declare commmand line options
my $seed;
my $db;
my $output_directory;
my $round_stop = 100;
my $vec_db;
my $expand = "B";
my $cap_param;
my $cross_param;
my $blast_param;
my $lce = 500;
my $qws = 0;
my $qwi = 0;
my $qpa = 10;
my $ml_chimera = 75;
my $out_files = "complete";
my $length_seed_stop = 10000000000; #just a big number;
my $put_last_contig_last_assembly = "yes";
my $minimum_length_block = 0;
my $minimum_length_block_id_percentage = 0;

my $help;

my $cont_html = 0;
my $number = 0;
my %contigs;
my %reads_blast;
my $im = "";



GetOptions( "s=s"      => \$seed,
	    "d=s"      => \$db,
	    "o=s"      => \$output_directory,
	    "r=s"      => \$round_stop,
	    "v=s"      => \$vec_db,
	    "e=s"      => \$expand,
	    "a=s"      => \$cap_param,
	    "c=s"      => \$cross_param,
	    "b=s"      => \$blast_param,
	    "l=s"      => \$lce,
	    "q=s"      => \$qws,
	    "p=s"      => \$qwi, 
	    "m=s"      => \$qpa,
	    "t=s"      => \$ml_chimera,
	    "x=s"      => \$out_files,
	    "f=s"      => \$length_seed_stop,
	    "g=s"      => \$put_last_contig_last_assembly,
	    "i=s"      => \$minimum_length_block,
	    "j=s"      => \$minimum_length_block_id_percentage,
	    "h"        => \$help
	    );


### print Help
my $help_print = "GenSeed version $version
\nUsage: genseed.pl -s <seed file> -d <database file> <optional parameters>\n
-h      Print this help

-s      Seed file name
-d      Database file name

Optional parameters:

-e      Expansion direction - WARNING - unidirectional expansion is only allowed for a single nucleic acid seed sequence
            l - expand only from the 5' end of the seed(s)
            r - expand only from the 3'end of the seed(s)
            b - expand to both ends of the seed(s) (default)
-f      Maximum length of the final consensus sequence (default: no limit)
-g      Drive final assembly with last-round consensus sequence (default: yes) 
        Warning: this option usually requires more processing time.
-i      Minimum length of the alignment block of the first BLAST search, in percentage relative of the length 
        of the seed sequence (default: 0)        
-j      Minimum identity percentage of the alignment block of the first BLAST search (default: 0)
-l      Length of contig end to be used as a seed in the next round (default: 500)
        (* and N characters are automatically removed from the sequence)	
-m      Maximum percentage of bad quality values in the sliding window (default: 10)
-o      Output directory (default: genseed_dir#)
-p      Minimum quality value for the sliding window (default: 0)
-q      Length of the contig quality verification window (default: 0)
-r      Maximum number of rounds (default: 100 rounds)
-t      Minimum length of misaligned ends to discard potentially chimeric sequences (default: 75)
-v      Vector file (vector sequence file to be used for masking purposes)
-x      Output files 
            complete - save the output of BLAST, CAP3 and generate HTML report files (default)
            simple - generate only HTML report and final assembly files

IMPORTANT: the parameters below MUST be defined under quotes. 
-a      CAP3 assembly parameters (example: \"-o 100 -p 95\")
-b      BLAST parameters (default: \"-F F -b 500 -e 1e-06\")
-c      Cross_match parameters for vector masking (default: \"-minmatch 12 -penalty -2 -minscore 20\")\n";


### Get command line options
if ($help){
    print "$help_print";
    exit;
}

die "$help_print\n" if !(($seed) && ($db));


if (($out_files !~ /complete/) && ($out_files !~ /simple/)){
    print "\nError\n-x valid options are complete or simple\n\n";
    exit;
}

if (!-e "$seed"){
    print "\nSeed file not found\n";
    exit;
}


if (!-e "$db"){
    print "\nDatabase file not found\n";
    exit;
}

if ($vec_db){
    if (!-e "$vec_db"){
	print "\n Vector database file not found\n";
	exit;
    }
}

if (($minimum_length_block < 0) || ($minimum_length_block > 100)){
    print "\nUsage:\n-i followed by a number between 0 and 100\n\n";
    exit;
}

if (($minimum_length_block_id_percentage < 0) || ($minimum_length_block_id_percentage > 100)){
    print "\nUsage:\n-j followed by a number between 0 and 100\n\n";
    exit;
}

### Verify seed type (DNA or protein)
my ($seed_type,$seed_quantity) = seed_verify($seed);
my $seed_orientation_control = 0;


if (($seed_quantity == 1 ) && ($seed_type==0)){
	###Compare the orientation of the former contig in regard to the seed sequence 
	$seed_orientation_control = 1
}

if (($seed_quantity > 1) && ($expand !~ /[Bb]/)){
    print "\n\nError! Unidirectional expansion is only allowed for a single nucleic seed sequence.\n\n";
    exit;
}

if (($seed_type == 1) && ($expand !~ /[Bb]/)){
    print "\n\nError! Unidirectional expansion is only allowed for a single NUCLEIC seed sequence.\n\n";
    exit;
}


###Save original seed type and name
my $original_seed_type = $seed_type;

my $seed_original_file_name = `basename $seed`;
chomp $seed_original_file_name;

### Verify if output directory has been defined
  # If not previously defined, a default directory with a increasing number is created
  #
my $contdir;
if (!$output_directory ){
    $contdir = 1;
    while (-d "genseed_dir$contdir"){
	$contdir ++;
    }
    $output_directory = "genseed_dir$contdir";
    system "mkdir $output_directory";
    system "mkdir $output_directory/CAP3_dir";
    system "mkdir $output_directory/image_dir";

    if ($out_files =~ m/complete/){
	system "mkdir $output_directory/BLAST_dir";
	system "mkdir $output_directory/fasta_dir";
    }
    
}
else{
    
    if (!(-d $output_directory)){
	system "mkdir $output_directory";
	system "mkdir $output_directory/CAP3_dir";
	system "mkdir $output_directory/image_dir";
	
	if ($out_files =~ m/complete/){
	    system "mkdir $output_directory/BLAST_dir";
	    system "mkdir $output_directory/fasta_dir";
	}
    }
    else {

	system "rm -rf $output_directory";
	system "mkdir $output_directory";
	system "mkdir $output_directory/CAP3_dir";
	system "mkdir $output_directory/image_dir";
	
	if ($out_files =~ m/complete/){
	    system "mkdir $output_directory/BLAST_dir";
	    system "mkdir $output_directory/fasta_dir";
	}
    }
}



system "cp $seed $output_directory";
chdir "$output_directory";

###Start genseed.log
my $genseed_log = "-s $seed -d $db ";
$genseed_log .= "-o $output_directory " if ($output_directory);
$genseed_log .= "-r $round_stop " if ($round_stop != 100);
$genseed_log .= "-v $vec_db " if ($vec_db);
$genseed_log .= "-e $expand " if ($expand !~ /B/);
$genseed_log .= "-a \\\"$cap_param\\\" " if ($cap_param);
$genseed_log .= "-c \\\"$cross_param\\\" " if ($cross_param);
$genseed_log .= "-b \\\"$blast_param\\\" " if ($blast_param);
$genseed_log .= "-l $lce " if ($lce != 500);
$genseed_log .= "-q $qws " if ($qws);
$genseed_log .= "-p $qwi " if ($qwi);
$genseed_log .= "-m $qpa " if ($qpa != 10);
$genseed_log .= "-t $ml_chimera " if ($ml_chimera != 75);
$genseed_log .= "-x $out_files " if ($out_files !~ /complete/);
$genseed_log .= "-f $length_seed_stop " if ($length_seed_stop != 10000000000);
$genseed_log .= "-g $put_last_contig_last_assembly " if ($put_last_contig_last_assembly !~ /yes/);
$genseed_log .= "-i $minimum_length_block " if ($minimum_length_block > 0);
$genseed_log .= "-j $minimum_length_block_id_percentage " if ($minimum_length_block_id_percentage > 0);
system "date >>genseed.log";
system "echo \"genseed.pl $genseed_log\" >>genseed.log";


$seed = `basename $seed`;
chomp $seed;

print "\nseed type: DNA\n" if ($seed_type==0);
system "echo \"seed type: DNA\" >>genseed.log" if ($seed_type==0);

print "\nseed type: Protein\n" if ($seed_type==1);
system "echo \"seed type: Protein\" >>genseed.log" if ($seed_type==1);


### Create headers for HTML files and verify if multifasta file needs to be formatted

new_germination ($db);

if (-e "../$db"){
   $db = "../$db";
}
elsif ($db =~ m/^..\//){
    $db = "../$db";
}

if ($vec_db){
    if (-e "../$vec_db"){
	$vec_db = "../$vec_db";
    }
    elsif ($vec_db =~ m/^..\//){
	$vec_db = "../$vec_db";
    }
}



### Get the fasta header(s) of seed file
my $header = get_header_seed ($seed);
my $header_original = $header;


##### Copy seed file to fasta_cap_final.fasta if seed type is DNA
if ($seed_type == 0){
    system "cat $seed >fasta_cap_final.fasta";
}


###for contig length stop
my $contig_length = "";

my $round = 0;
my $one_read_flag = 0;

### Main loop
MAIN_LOOP:
while (1){

    $round++;


    ## finish by round 
    if ($round > $round_stop){
	finishing($db, $cap_param,$header_original,$round);
    }
    
    print "\n####   Round $round   \####\n";
    system "echo \"####   Round $round   ####\" >>genseed.log";

    #### Blast seed sequence file against db###
    if ($seed_type == 0){
	blast ($seed,$db,$blast_param,"blastn");
	system "cp blast_output BLAST_dir/blast_$round.out"  if ($out_files =~ /complete/);
    }
    else{
	#### Blast protein seed file against db
	blast ($seed,$db,$blast_param,"tblastn");
	system "cp blast_output BLAST_dir/blast_$round.out"  if ($out_files =~ /complete/);
    }
    
    #### Parse blast output file to extract read information
    my $control_end_blast = parser_blast($db,$seed,$round);

       
    ### If no new new sequence is found, start the finishing process
    if ($control_end_blast == 1){

	if ($round > 1){
	    print "No new read!!!\n";
	    system "echo \"No new read!!!\" >>genseed.log";
	    finishing($db, $cap_param,$header_original,$round);
	}
	else{
	    print "No match with selected database!\nPlease choose another seed.\n\n";
	    system "echo \"No match with selected database!\nPlease choose another seed\" >>genseed.log";
	    exit;
	}
    }


    ############ IF only one read is selected #####################
    #if (($round == 1) && ($seed_type == 1)){
    if ($round == 1){
	my $number_of_selected_reads = `grep \">\" accepted_reads.fasta|wc -l`;
	
	if ($number_of_selected_reads == 1){
	    
	    $one_read_flag = 1;
	    system "cp accepted_reads.fasta one_read.fasta";
	    
	    print "\nOnly one read has been selected!\nUsing this read for Blast similarity search\n\n";
	    
	    system "cp accepted_reads.fasta accepted_reads.temp";
	    system "mv accepted_reads.fasta fasta_cap.fasta";
	    
	    my $seed_name_temp = $seed;
	    $seed = "accepted_reads.temp";
	    
	    blast ($seed,$db,$blast_param,"blastn");
	    system "cp blast_output BLAST_dir/blast_a_$round.out"  if ($out_files =~ /complete/);
	    
	    my $control_end_blast = parser_blast($db,$seed,$round);
	    
	    $seed = $seed_name_temp;
	    system "cat accepted_reads.temp >>cap_input";
	    system "rm accepted_reads.temp";

	    my $number_of_selected_reads = `grep \">\" accepted_reads.fasta|wc -l`;

	    if ($number_of_selected_reads == 1){		
		print "Only one DNA read matches your seed sequence, but it does not present 
		any overlapping reads!\n Final contig corresponds to this sequence\n\n";
		system "echo \"Only one DNA read matches your seed sequence, but it does not present 
		any overlapping reads! Final contig corresponds to this sequence.\" >>genseed.log";

		system "mv one_read.fasta final_contigs.fasta";
		exit;
	    }
	}
    }
    ############  END IF only one read is selected #####################


    ##### Remove potential chimeric reads
    my $round_after = $round-1;
    remove_chimeric_reads ("accepted_reads.fasta","contig_$round_after.fasta",$ml_chimera,$blast_param) if ($round > 1);


    #### Mask vector bases with Cross_match
    cross_match("accepted_reads.fasta",$vec_db,$cross_param,"0") if ($vec_db =~ m/./);
    
    #### Run CAP3 on the new set of read(s) plus the last contig sequence
    cap3($seed,$cap_param,$round,$seed_type);
    
    system "cp cap_input.cap.ace CAP3_dir/cap_$round.ace";

    if ($out_files =~ /complete/){
	system "cp cap_input.cap.contigs CAP3_dir/cap_$round.contigs";
	system "cp cap_input.cap.singlets CAP3_dir/cap_$round.singlets";
	system "cp cap_input CAP3_dir/cap_$round.input";
    }

    #### Parse CAP3 output file and compare current and former contigs (length and read composition)
    my $control_end;
    my @fastas;
    my $contig_number;
    my $rep_seed_position;

    ($control_end,$contig_number,$rep_seed_position,$header,$contig_length,@fastas) = parser_cap3($db, $header,"cap_input.cap.ace", "0",$expand,$contig_length,$round,$seed,$seed_orientation_control,$seed_type);
    
    system "cp contig_round.fasta contig_$round.fasta";
    if ($out_files =~ /complete/){
	system "cp fasta_cap.fasta fasta_dir/seed_$round.fasta";
    }
    
    #### If the current contig presents the same number of bases or reads, stop the main loop 
    #### and start the finishing process 
    #### If contigs are different, send contig assembly information to the report sub-routine
    

    for my $aux (@$contig_length){
	if ($aux > $length_seed_stop){
	    print "A contig reached the maximum number of bases!\n";
	    $control_end = 0;
	}

    }

    if ($control_end == 0){
	report($round,"0",$contig_number,$seed,\@fastas,$seed_type,$rep_seed_position);
	finishing($db, $cap_param,$header_original,$round);
    }
    else{
	report($round,"0",$contig_number,$seed,\@fastas,$seed_type,$rep_seed_position);
    }

    system "rm contig_round.fasta";
    system "rm temporary_seed_fasta.fasta*" if (-e "temporary_seed_fasta.fasta");

    $seed = "fasta_cap.fasta";

    ##### Seed type is now a DNA composed by the ends of the consensus sequence
    $seed_type = 0;
}

exit;

####################
### Sub-routines ###
####################



### Open seed file to get the fasta headers of the seed(s)
sub get_header_seed{

    my $seed = shift;

    my @header = `grep ">" $seed`;
    my @header1;

    for (@header){
	$_ =~ m/>(\S+)/;
	push @header1, $1;
    }

    ##### Open report.html to print seed name(s)
    open (HTML, ">>report.html");

    my $seed_number = @header1+0;
    print HTML "Number of seeds: $seed_number<BR>\n";

    print HTML "Seed name(s):<BR>";
    
    for (@header1){
	print HTML "- $_<BR>\n";
    }
    
    close(HTML);
    return \@header1;
}

sub new_germination {

    my $db = shift;

    #####create report.html
    open (HTML, ">report.html");
    print HTML "<html>
<head>
<title>GenSeed - Graphic report<\/title>
<\/head>
<body>

<p>GenSeed - Graphic report<\/p>
";
    close(HTML);

    #####create list_of_reads.html
    open (HTML_READS, ">list_of_reads.html");
    print HTML_READS "<html>

<head>
<title>GenSeed - Read report<\/title>
<\/head>

<body>

<p>GenSeed - Read report<\/p>

The lists below present all sequence reads that were included in the assembly of <BR>
each round.<BR><BR>

";
    close(HTML_READS);

    ##### Create seed_contigs.html
    open (HTML_CONTIGS, ">seed_contigs.html");
    print HTML_CONTIGS "<html>

<head>
<title>GenSeed - Contig reporter<\/title>
<\/head>

<body>

<p>GenSeed - Contig report<\/p>

The sequences below represent the different contigs obtained in each one of 
the assembly rounds. These contigs necessarily contain the seed(s) sequence(s)
when the seed was a nucleotide sequence. For protein seeds, the seed 
sequence is only listed from the second round.
<BR>
";
    close(HTML_CONTIGS);

    chdir "../";

    ### Verify if formatdb has already processed the multifasta file and created the corresponding index files
    ### If index files are present, formatdb will not run again
    if (-e ("$db.nhr") && -e ("$db.nin") && -e ("$db.nsd") && -e ("$db.nsi") && -e ("$db.nsq")){ 
    }
    elsif (-e ("$db.nal")){

	open (DB_DATA, "$db.nal");
	my @files_db;
	while (my $formatdb_nal = <DB_DATA>){
	    if ($formatdb_nal =~ m/^DBLIST (.*)/){
		@files_db = split (/ /,$1);
		last;
	    }
	}

	close (DB_DATA);

	foreach my $aux (@files_db){
	    if (-e ("$aux.nhr") && -e ("$aux.nin") && -e ("$aux.nsd") && -e ("$aux.nsi") && -e ("$aux.nsq")){ 
		
	    }
	    else{
		system "formatdb -o T -p F -i $db 2>>genseed.log";
		system "rm formatdb.log";
	    }
	}
    }
    else{
	system "formatdb -o T -p F -i $db 2>>genseed.log";
	system "rm formatdb.log";
    }

    chdir "$output_directory";

    ### Eliminate old files
    system "rm fasta_cap.fasta" if (-e "fasta_cap.fasta");
    system "rm fasta_cap_final.fasta" if (-e "fasta_cap_final.fasta");

    return();
}

### Run blast 
sub blast {

    my $file = shift;
    my $blastdb = shift;
    my $blast_par = shift;
    my $blast_program = shift;

    ##### Verify if Blast parameters were not defined in command line
    if (!$blast_par){
	$blast_par = "-F F -b 500 -e 1e-06";
    }

    system "blastall -p $blast_program -i $file -d $blastdb -o blast_output -m8 $blast_par 2>>genseed.log";
    return;
}



### Parse Blast output file

sub parser_blast{

    my $blastdb = shift;
    my $seed_file_name = shift;
    my $round = shift;

    my %fastas_header;
    

    if (($round == 1) && (($minimum_length_block>0) || ($minimum_length_block_id_percentage>0))){
	open (BLAST, "$seed_file_name") or die "Cannot open $seed_file_name!";
	
	my $cont = 0;
	my $fasta_header;
	while (<BLAST>){
	    chomp;
	    if ($_ !~ /^>(\w+)/){
		 chomp $_;
		 my $aux = $_;
		 $aux =~ s/\s//g;
		 $cont += length $aux;
	    }
	    else{

		if ($_ =~ /^>.*/){
		    $fastas_header{$fasta_header} = $cont;
		    $_ =~ m/^>(\w+)/;
		    $fasta_header = $1;
		    $cont = 0;
		}
	    }
	}
	
	$fastas_header{$fasta_header} = $cont;
    }
    my $control_end = 0;

    ##### Open Blast output file
    open (BLAST, "blast_output") or die "Cannot open blast_output!";

    ##### Remove temporary files
    system "rm accepted_reads.fasta*\n" if (-e "accepted_reads.fasta");

    while (<BLAST>){
	
	$_ =~ m/(\S+)\s+(\S+)\s+([\d\.]+)\s+(\d+)\s+\d+\s+\d+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+/;
	my $query_name = $1;
	my $name_subject = "\"$2\"";
	my $name_subject1 = $2;

	my $identity = $3;
	my $length_block = $4;

	if (($round == 1) && (($minimum_length_block>0) || ($minimum_length_block_id_percentage>0))){

	    if ($length_block < (($minimum_length_block/100)*$fastas_header{$query_name})){
		next;
	    }
	    else{

		if ($identity < $minimum_length_block_id_percentage){
		    next;
		}
		
	    }
	}

	### Avoid read duplication 
	if (exists $reads_blast{$name_subject1}){
	    next;
	}
	else{
	    $reads_blast{$name_subject1} = 1;
	    $control_end = 1;
	}

	system "fastacmd -d $blastdb -s $name_subject >>accepted_reads.fasta";
    }

    ### If no new read is found, return a 1 and go to the finishing process
    return (1) if ($control_end == 0);
    system "cp accepted_reads.fasta cap_input";
    return;
}

sub remove_chimeric_reads {

    my $file = shift;
    my $contig_file = shift;
    my $quant_elemination_bases = shift;
    my $blast_par = shift;

    my $reads_size = reads_size ($file);
    my %read_size = %{$reads_size};    

    my $sequence_size = reads_size ($contig_file);
    my %sequence_size = %{$sequence_size};    


    my %read_alg_size;


    system "formatdb -o T -p F -i $file 2>>genseed.log";

    if (!$blast_par){
	$blast_par = "-F F -b 500 -e 1e-06";
    }

    my @blast = `blastall -p blastn -i $contig_file -d $file $blast_par -m8 2>>genseed.log`; 

    for my $aux (@blast){

	my @blast_aux = split (/\s+/,$aux);

	my $control_start_contig = 0;
	if ($blast_aux[6] > $quant_elemination_bases){
	    $control_start_contig = 1;
	}

	my $control_end_contig = 0;
	if ($blast_aux[7] < $sequence_size{$blast_aux[0]} - $quant_elemination_bases){
	    $control_end_contig = 1;
	}
	
	if ($blast_aux[9] < $blast_aux[8]){
	    if ($read_alg_size{$blast_aux[1]} < $blast_aux[3]){
		$read_alg_size{$blast_aux[1]} = $blast_aux[3];

		#####Verify if the end of the read 
 		if  (($control_start_contig == 1) && ($control_end_contig == 1) && ($blast_aux[8] < $read_size{$blast_aux[1]} - $quant_elemination_bases)){
		    next;
		}
		elsif  (($control_start_contig == 1) && ($control_end_contig == 1) && ($blast_aux[9] > $quant_elemination_bases)){
		    next;
		}
		else{
		    system "fastacmd -d $file -s \"$blast_aux[1]\" >>temporary_fasta.fasta 2>>genseed.log";
		}
	    }
	}
	else{

	    if ($read_alg_size{$blast_aux[1]} < $blast_aux[3]){
		
		$read_alg_size{$blast_aux[1]} = $blast_aux[3];

		if  (($control_start_contig == 1) && ($control_end_contig == 1) && ($blast_aux[9] < $read_size{$blast_aux[1]} - $quant_elemination_bases)){
		    next;
		}
		elsif  (($control_start_contig == 1) && ($control_end_contig == 1) && ($blast_aux[8] > $quant_elemination_bases)){
		    next;
		}
		else{
		    system "fastacmd -d $file -s \"$blast_aux[1]\" >>temporary_fasta.fasta 2>>genseed.log";
		}
	    }
	}
    }

    system "rm $file.n*";
    system "cp temporary_fasta.fasta $file" if (-e "temporary_fasta.fasta");
    system "cp temporary_fasta.fasta fasta_cap.fasta" if (-e "temporary_fasta.fasta");
    system "mv temporary_fasta.fasta cap_input" if (-e "temporary_fasta.fasta");

    return ();
}

sub reads_size {

    my $file = shift;

    my $fasta_name;
    my %fasta_length;


    open (DATA, "$file");
    
    while (my $line = <DATA>){
	if ($line =~ /^>lcl\|(\S+)/){
	    $fasta_name = $1;
	}
	elsif ($line =~ /^>(\S+)/){
	    $fasta_name = $1;
	}
	else{
	    chomp $line;
	    $fasta_length{$fasta_name} += length $line;
	}
    }
    close (DATA);
    return (\%fasta_length);
}



### If a vector file has been defined, run Cross_match to mask vector bases
sub cross_match {

    my $file = shift;
    my $vec_base = shift;
    my $cross_par = shift;
    my $final_control = shift;

    ##### Verify if Cross_match parameters were not defined in command line
    if (!$cross_par){
	$cross_par = "-minmatch 12 -penalty -2 -minscore 20";
    }

    system "cross_match $cross_par -screen $file $vec_base >cross_temp_genomic 2>>/dev/null";

    system "cp accepted_reads.fasta.screen cap_input" if ($final_control == 0);
    system "cp fasta_cap_final.fasta.screen fasta_cap_final.fasta" if ($final_control == 1);

    return;
}


### Run CAP3 program

sub cap3{

    my $file = shift;
    my $cap_par = shift;
    my $round = shift;
    my $seed_type = shift;
    $round--;

    `cat $file >> cap_input`if (($round == 0) && ($seed_type == 0));
    `cat contig_$round.fasta >> cap_input` if ($round > 0);

    my $cap_reads = `grep ">" cap_input|wc -l`;
    $cap_reads = $cap_reads-$seed_quantity if ($seed_type == 0);
    chomp $cap_reads;
    print "Total # of reads for CAP3: $cap_reads\n";
    system " echo \"Total # of reads for CAP3: $cap_reads\" >>genseed.log";

    system "cap3 cap_input $cap_par >cap3_temp_genomic 2>>genseed.log";
    
    return;
}


### Parse CAP3 files and verify if the seed-containing contig(s) presents newly incorporated bases and/or reads
### Store information for report files

sub parser_cap3 {
    
    my $blastdb = shift;
    my $header = shift;
    my $ace_file = shift;
    my $last_run = shift;
    my $expand = shift;
    my $contig_length = shift;
    my $round = shift;
    my $seed_file_name = shift;
    my $seed_orientation_control = shift;
    my $seed_type = shift;

    my $control_finish = 0;
    my @fastas = (); 
    
    open (DATA, "$ace_file");

    my %sequence;
    my %quality;
    my %contig;
    my %contigs_seed;
    my %rep_seed_position;
    my $contig_number = "";

    while (<DATA>){

	#####extract contig number and contig senquence
	if ($_ =~ m/CO\s+Contig(\d+)/){
	    $contig_number = $1;
	    my $aux_sequence = $_;
	    my $sequence = "";

	    while ($aux_sequence =~ m/\w+/){
		$aux_sequence = <DATA>;
		chomp $aux_sequence;
		$aux_sequence =~ s/\s//g;
		$sequence .= $aux_sequence;
	    }
	    $sequence{$contig_number} = $sequence;
	}
	
	##### Extract contig base quality
	if ($_ =~ m/BQ/){

	    my $aux_quality = $_;
	    my $quality = "";

	    while ($aux_quality =~ m/\w+/){
		$aux_quality = <DATA>;
		chomp $aux_quality;
		$quality .= $aux_quality;
	    }

	    $quality{$contig_number} = $quality;
	}

	##### Extract contig reads
	if ($_ =~ m/AF\s+(.*)\s+\w\s+[-\d]+/){

	    my $fasta_name = $1;
	 
	    if ($fasta_name !~ m/contig_\d+_seed_program/){
		$contig{$contig_number}{$fasta_name} = 1;
	    }

	    ##### Verify if seed DNA sequence is present in one of the reads
	    if ($seed_type == 0){
		foreach my $aux (@$header){
		    
		    my $fasta_name_aux = $fasta_name;
		    $fasta_name_aux =~ s/\|/\./g;
		    $aux =~ s/\|/\./g;
		 
		    if ($fasta_name_aux =~ /$aux/){
			delete  $contig{$contig_number}{$fasta_name};
			
			if (($round > 1) && ($last_run == 0)) {
			    if ($aux =~ m/contig_\d+_seed_program_(.*)/){
				next if ($contigs_seed{$contig_number} =~ /$1/);
				$contigs_seed{$contig_number} .= "_$1";		    
			    }
			    else{
				next if ($contigs_seed{$contig_number} =~ /$aux/);
				$contigs_seed{$contig_number} .= "_$aux";		    
			    }
			}
		    }
		}
	    }
	}
    }
    close (DATA);

    #### For first or last Round
    if (($round == 1) || ($last_run == 1)){

	#### Create a temporary CAP3 contigs database file
	open (FASTA_TEMP, ">temporary_seed_fasta.fasta");

	for my $key (keys %sequence){
	    print FASTA_TEMP ">$key\n$sequence{$key}\n";
	}
	close (FASTA_TEMP);

	system "formatdb -o T -p F -i temporary_seed_fasta.fasta 2>>genseed.log\n";

	my @blast_out = ();

	if ($seed_type == 1){
	    @blast_out = `blastall -p tblastn -d temporary_seed_fasta.fasta -i $seed_original_file_name -m8 -e 10000 2>>genseed.log`;
	}
	
	if ($seed_type == 0){
	    @blast_out = `blastall -p blastn -d temporary_seed_fasta.fasta -i $seed_original_file_name -m8 -e 10000 2>>genseed.log`;
	}

	system "rm temporary_seed_fasta.fasta*" if (-e "temporary_seed_fasta.fasta");
	
	
	my $s_id_old;
	my $q_id_old;
	my $contig_seed;
	my $contig_seed_name;
	my $percentage_old = 0;
	my $alignment_length_old = 0;
	my $control_first = 0;
	my $rep_seed_position;
	
	my $q_id;
	my $q_begin;
	my $q_end;
	my $s_begin;
	my $s_end;

	


	for my $blast_out_aux (@blast_out){

	    #######print stderr "$blast_out_aux";
	    if ($blast_out_aux =~ m/No hits found/){
		last;
	    }

	    if ($blast_out_aux =~ m/(\S+)\s+(\S+)\s+([\d\.]+)\s+(\d+)\s+\d+\s+\d+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+/){

		$q_id = $1;
		my $s_id = $2;
		my $percentage = $3;
		my $alignment_length = $4;

		$q_begin = $5;
		$q_end = $6;
		$s_begin = $7;
		$s_end = $8;
		

		if ($control_first == 0){
		    $q_id_old = $q_id;
		    $control_first = 1;
		}

		if ($q_id eq $q_id_old){
		    if ($alignment_length > $alignment_length_old){
			if ($percentage > $percentage_old){

			    $percentage_old = $percentage;
			    $alignment_length_old = $alignment_length;
			    
			    $contig_seed = $s_id;
			    $contig_seed_name = $q_id;
			    if ($s_begin < $s_end){
				$rep_seed_position = "$s_begin $s_end $q_id##";
			    }
			    else{
				$rep_seed_position = "$s_end $s_begin $q_id##";
			    }
			}
		    }
		}
		else{
		    
		    $contigs_seed{$contig_seed} .= "_$contig_seed_name";

		    $rep_seed_position{$contig_seed} .= $rep_seed_position;
    
		    $percentage_old = $percentage;
		    $alignment_length_old = $alignment_length;
		    
		    $contig_seed = $s_id;
		    $contig_seed_name = $q_id;
		    if ($s_begin < $s_end){
			$rep_seed_position = "$s_begin $s_end $q_id##";
		    }
		    else{
			$rep_seed_position = "$s_end $s_begin $q_id##";
		    }

		}

		$q_id_old = $q_id;
	    }
	}

	if ($contig_seed){
	    $contigs_seed{$contig_seed} .= "_$contig_seed_name";
	    $rep_seed_position{$contig_seed} .= $rep_seed_position;
	}
    }


    ##### Verify if seed senquence is found in any contig
    my $seed_contig_control = 0;
    for (keys %contigs_seed){

	#print "DEBUG:1:# $_ #$contigs_seed{$_}  # $rep_seed_position{$_}\n";
	
	$seed_contig_control = 1;
    }


    ##### If seed senquence is not present in any contig, then start the finishing process or 
    ##### quit if this is the first round
    if (($seed_contig_control == 0) && ($last_run == 0)){	
	
	if (%contigs > 0){
	    print "Fasta seed is not present in any contig!\nFinishing...\n\n";
	    system "echo \"Fasta seed is not present in any contig!\nFinishing...\" >>genseed.log";
	    
	    finishing($db,$cap_param,$header_original,$round);
	}
	else{
	    print "Fasta seed is not present in any contig!\nExit\n\n";
	    system "echo \"Fasta seed is not present in any contig!\nExit\" >>genseed.log";

	    if (($round==1) && ($one_read_flag == 1)){
		system "mv one_read.fasta final_contigs.fasta";
		if (-e "blast_output"){
		    system "rm blast_output";
		}
		if (-e "fasta_cap.fasta"){
		    system "rm fasta_cap.fasta";
		}
		if (-e "fasta_cap_final.fasta"){
		    system "rm fasta_cap_final.fasta";
		}

	    }
	    
	    system "rm cap_input* cap3_temp_genomic $seed_original_file_name accepted_reads.fasta";
	    system "rm cross_temp_genomic" if (-e "cross_temp_genomic");
	    
	    exit;
	}
    }


    ##### IF LAST ROUND, return seed-containing contig(s) and contig reads
    if ($last_run == 1){
	for (keys %contigs_seed){
	    
	    my $sequence = $sequence{$_};
	    $sequence =~ s/\*//g;
	    $sequence =~ s/\s//g;
	    $sequence =~ s/\n//g;

	    my $length_seq_contig_aux = length $sequence;

	    my $seed_name_aux = $contigs_seed{$_};
	    $seed_name_aux =~ s/_//;
	    $seed_name_aux =~ s/\|/./g;
	    
	    my $name_aux = $seed_name_aux;
	    $name_aux =~ s/\|/./g;
	    print "Length of the final seed-positive contig $name_aux: $length_seq_contig_aux\n";
	    system "echo \"Length of the final seed-positive contig $name_aux: $length_seq_contig_aux\" >>genseed.log";
	    
	}


	return (\%contigs_seed, \%contig, \%rep_seed_position, $seed_contig_control);
    }

    ##### If not last round, extract contig ends for the next round and store new reads 

    open (FASTA_OUTPUT, ">fasta_cap.fasta");
    
    ##### Open files temporary Genseed file
    open (CONTIG_OUTPUT, ">contig_round.fasta");

    ##### Open files containing formatted versions of the consensus sequences
    open (CONSENSUS_OUTPUT, ">fasta_dir/consensus_$round.fasta"); 

    
    open (HTML_READS, ">>list_of_reads.html");
    open (HTML_CONTIGS, ">>seed_contigs.html");


    print HTML_READS "<b>Round #$round<\/b><BR><BR>\n";
    print HTML_CONTIGS "<b>Round #$round<\/b><BR><BR>\n";

    #print CONSENSUS_OUTPUT "Round #$round\n\n";



    my @contig_number = ();
    my @header = ();
    my @contig_length = ();
    my $control_length_contigs_final = 0;

    ##### Separate new reads
    for my $aux (keys %contigs_seed){

	$sequence{$aux} = contig_orientation ($sequence{$aux}, $seed_original_file_name) if ($seed_orientation_control == 1);

	push @contig_number, "$aux\_$contigs_seed{$aux}";

	for (keys %{$contig{$aux}}){
	    next if ($_ =~ m/contig_\d+_\d+_\d+/);
	    
	    my $header_string = $_;

	    ##### If the read exists...
	    if (exists $contigs{$header_string}){
		#not new!!! 
		push @fastas, $header_string;
	    }
	    ##### If the read is new, add it to the final read list
	    else{
		# It is new...
		print HTML_READS "$header_string<BR>\n";
		push @fastas, $header_string;
		$contigs{$header_string} = 1;
		$control_finish = 1;
	    }
	}

	##### Create base quality and sequence arrays for quality inspection and
	##### report.html file
	my @quality = split (/\s/, $quality{$aux});
	my @sequence_aux = split (//, $sequence{$aux});
	my $seq_contig_aux = $sequence{$aux};
	$seq_contig_aux =~ s/\*//g;

	my $seq_html = $seq_contig_aux;

	my $seed_name_aux = $contigs_seed{$aux};
	$seed_name_aux =~ s/_//;
	$seed_name_aux =~ s/\|/./g;

	print HTML_CONTIGS "<FONT FACE=\"Courier\">>$seed_name_aux - extended sequence - Contig$aux<BR>\n";
	print CONSENSUS_OUTPUT ">$seed_name_aux - extended sequence - Contig$aux\n";

	### Temporary GenSeed file
	print CONTIG_OUTPUT ">contig_$aux\_seed_program$contigs_seed{$aux}\n";

	##### Format contig sequences to 60 bp per line
	while ($seq_html){
	    my $aux = substr ($seq_html,0,60,"");
	    print HTML_CONTIGS "$aux<BR>\n";	    
	    print CONSENSUS_OUTPUT "$aux\n";
	    print CONTIG_OUTPUT "$aux\n";	    		
	}

	print HTML_CONTIGS "</FONT>\n";

	my $name_aux = $seed_name_aux;
	$name_aux =~ s/\|/./g;


	my $seq_contig_aux_length = $seq_contig_aux;
	
	$seq_contig_aux_length =~ s/\*//g;
	$seq_contig_aux_length =~ s/\s//g;
	$seq_contig_aux_length =~ s/\n//g;

	my $length_seq_contig_aux = length $seq_contig_aux_length;

	print "Length of the seed-contig $name_aux: $length_seq_contig_aux\n";
	system "echo \"Length of the seed-contig $name_aux: $length_seq_contig_aux\" >>genseed.log";
	
	my $control_length_contigs = 0;
	
	##### Verify contig length
	if ($round > 1 ){
	    for my $aux_length (@$contig_length){
		if (length $seq_contig_aux == $aux_length){
		    $control_length_contigs = 1;
		}
	    }
	}
	
	if ($control_length_contigs == 0){
	    $control_length_contigs_final = 1;
	}

	push @contig_length, length $seq_contig_aux;
	push @header, "contig_$aux\_seed_program$contigs_seed{$aux}";

	##### Prepare contig sequence to sliding window quality inspector
	my $temporary_counter = 0;
	my @quality_temporary_array = ();
	for (@sequence_aux){

	    if ($_ =~ m/\*/){
		push @quality_temporary_array,$qwi;
		next;
	    }
	    push @quality_temporary_array, $quality[$temporary_counter];
	    $temporary_counter++;
	}

	my $tam_quality_temporary_array = @quality_temporary_array+0;

	my $initial_pos;

	##### Quality windows inspector (qwi)
	    # This routine uses a sliding window to verify if the sequence ends present stretches with base 
	    # qualities above a user-defined value. Since FASTA sequences are used, CAP3 ascribes a quality
	    # value of 20 to all bases. However, according to the read redundancy, CAP3 gives quality values 
	    # to the consensus sequence that can be as high as 90.
	    #
	    # Qwi is used to discard end bases with low consensus quality (usually a value between 20 and 90)
	    # and to define the coordinates (starting from the 5' end) from which the sequence ends 
	    # will be extracted and used for the next round.
	
	##### Define starting coordinate of 5' end
      quality_loop:
	for (my $a = 0; $a <= $tam_quality_temporary_array; $a++){
	    my $end_b_loop = $a+$qws;
	    my $bad_quality_value = 0;

	    for (my $b = $a; $b <= $end_b_loop; $b++){
		if ($quality_temporary_array[$b] < $qwi){
		    
		    $bad_quality_value++;
		    
		    if (($bad_quality_value/$qws) > ($qpa/100)){
			next quality_loop;
		    }
		}
	    }
	    $initial_pos = $a;
	    
	    last;
	}
	
	##### Create 5' end seed	
	my $ini_contig  = substr ($sequence{$aux}, $initial_pos, $lce); 
	$ini_contig =~ s/\*|N|n//g;
	my $initial_pos_500 = $initial_pos + $lce;
	my $final_pos;

	##### Define starting coordinate of 3' end
      quality_loop_end:
	for (my $a = $tam_quality_temporary_array; $a >= $initial_pos_500; $a--){

	    my $tam_quality_aux = $a - $qws;	
	    my $bad_quality_value = 0;

	    for (my $b = $a; $b >= $tam_quality_aux; $b--){
		if ($quality_temporary_array[$b] < $qwi){

		    $bad_quality_value++;
		    
		    if (($bad_quality_value/$qws) > ($qpa/100)){
			next quality_loop;
		    }
		}
	    }
	    $final_pos = $a-$lce;
	    last;
	}
	
	##### Create 3' end seed
	my $end_contig  = substr ($sequence{$aux}, $final_pos, $lce); 
	$end_contig =~ s/\*|N|n//g;

	my $final_pos_500 = $final_pos + $lce;

	if ($expand =~ m/[BbLl]/){
	    print FASTA_OUTPUT ">contig_$aux\_$initial_pos\_$initial_pos_500\_start
$ini_contig
";
	}
	
	if ($expand =~ m/[BbRr]/){
	    print FASTA_OUTPUT ">contig_$aux\_$final_pos\_$final_pos_500\_end
$end_contig
";
	}
    }
    close (FASTA_OUTPUT);
    close (CONTIG_OUTPUT);

    
    print HTML_READS "<BR>\n";
    print HTML_CONTIGS "<BR>\n";
    print CONSENSUS_OUTPUT "\n";

    close (HTML_CONTIGS);
    close (HTML_READS);
    close (HTML_READS);

    my $quant_read_fasta_final = 0;

    ##### Quantify reads for the last round
    for (keys %contigs){
	$quant_read_fasta_final++;
    }
    
    print "Accumulative number of reads: $quant_read_fasta_final\n";
    system "echo \"Accumulative number of reads: $quant_read_fasta_final\" >>genseed.log";

	
    ##### Finish walking cycling by number of reads
    if ($control_finish == 0){
	print "End! Same number of reads\n";
	system "echo \"End! Same number of reads\" >>genseed.log";
    }

    ##### Finish walking cycling by number of bases
    if ($control_length_contigs_final == 1){
	$control_finish = 1;
    }
    else{
	$control_finish = 0;
	print "The end! Contigs have the same size\n";
	system "echo \"The end! Contigs have the same size\" >>genseed.log";
    }

      return ($control_finish, \@contig_number,\%rep_seed_position,\@header,\@contig_length, @fastas);
}


sub contig_orientation {
    my $sequence = shift;
    my $seed_file_name = shift;
    my $blast_seed;

    $sequence =~ s/\*//g;

    open (CONTIG_FILE, ">aux_contig_file.fasta") or die "Unable to open aux_contig_file.fasta" ;
    print CONTIG_FILE ">contig_sequence\n$sequence\n";
    close (CONTIG_FILE);

    my @blast = `bl2seq -p blastn -i $seed_file_name -j aux_contig_file.fasta -e 1e-06 -F F -D 1 2>>genseed.log`;

    for my $aux (@blast){

	my @blast_aux = split (/\s+/,$aux);

	if ($blast_aux[9] < $blast_aux[8]){
	    if ($blast_aux[3] > $blast_seed){
		    $blast_seed = $blast_aux[3];
	    }
	}
	else{
	    if ($blast_aux[3] > $blast_seed){
		$blast_seed = "n";
	    }
	}
    }

    system "rm aux_contig_file.fasta";
    
    if ($blast_seed !~ /n/){
	    my $aux_sequence = reverse $sequence;
	    $aux_sequence =~ tr/AaTtGgCc/TtAaCcGg/;
	    return ($aux_sequence);

    }
    else{
	return ($sequence);
    }
}


##### Create HTML report files
sub report{

    my $round = shift;
    my $final = shift;
    my $contig_number = shift;
    my $seed = shift;
    my $fastas = shift;
    my $seed_type = shift;
    my $rep_seed_position = shift;

    my @contig_aux = ();
    my $map;
    my $head;
    my $graphic_contig_number;

    open (HTML, ">>report.html");

    if ($final == 0){

	print HTML "<p> Round #$round <\/p>\n<ol>\n";
	
	my $quant_reads = @$fastas+0;

	print HTML "</ol>\n";
	
        #################### Create graphical outputs for intermediate steps
	($map,$head,$graphic_contig_number) = gif_generation ("CAP3_dir/cap_$round.ace",$contig_number,$round,$seed,"",$seed_type,$rep_seed_position);

	my @ls_graphical_representation = split (/\#\#/,$graphic_contig_number);
	my $cont_map = 0;
	for my $aux_graphical_representation (@ls_graphical_representation){

	    $aux_graphical_representation =~ m/(\d+)/;

	    system "mv image_dir/graphical_representation$1.png image_dir/graphic_$round\_$1.png";

	    my $final_graphic_number = $1;
	    
	    $cont_html++;

	    print HTML "
<head>
@$head[$cont_map]
<\/head>
<map name=\"mapname$cont_html\">
@$map[$cont_map]
<\/map>
<img src=\"image_dir/graphic_$round\_$final_graphic_number.png\" USEMAP=\"#mapname$cont_html\" border=1><BR><BR>\n\n";
	    $cont_map++;
	}
        ########################## END graphics


	print HTML "Total number of reads incorporated in this round: $quant_reads<BR>\n";
	
	my $quant_read_fasta_final = 0;
	
	for (keys %contigs){
	    $quant_read_fasta_final++;
	}
	
	print HTML "Accumulative number of reads: $quant_read_fasta_final<BR>\n";
	
    }
    else{

	print HTML "<p> Last round<\/p>\n<ol>\n";
	print HTML "</ol>\n";

	my %aux_reads = %$fastas;

	open (HTML_READS, ">>list_of_reads.html");
	print HTML_READS "### FINAL CONTIG(S) READS ###<BR>\n";

	my $quant_reads = 0;

	my %contig_number_aux = %{$contig_number};

	#for my $aux (keys %{$contig_number}){
	for my $aux (keys %contig_number_aux){

	    for my $aux_html (keys %{$aux_reads{$aux}}){
		print HTML_READS "<li>$aux_html<\/li>\n";
		$quant_reads++;
		push @contig_aux, "$aux\__$contig_number_aux{$aux}";
	    }
	}
    	
	print HTML_READS "</ol>\n";
	close (HTML_READS);


        #################### Create graphical outputs for the final assembly step
	if ($round =~ m/\d+/){
	    ($map,$head,$graphic_contig_number) = gif_generation ("CAP3_dir/cap_$round.ace",$contig_number,$round,$seed,"",$seed_type,$rep_seed_position);
	}
	else{
	    ($map,$head,$graphic_contig_number) = gif_generation ("fasta_cap_final.fasta.cap.ace",\@contig_aux,1,$seed,$final,$seed_type,$rep_seed_position);
	}

	my @ls_graphical_representation = `ls  image_dir/graphical_representation*.png`;

	chomp @ls_graphical_representation;
	
	my $cont_map = 0;
	for my $aux_graphical_representation (@ls_graphical_representation){
	    $aux_graphical_representation =~ m/graphical_representation(\d+).png/;
	    system "mv $aux_graphical_representation image_dir/graphic_$round\_$1.png";
	    
	    $cont_html++;
	    print HTML "
<head>
@$head[$cont_map]
<\/head>
<map name=\"mapname$cont_html\">
@$map[$cont_map]
<\/map>
<img src=\"image_dir/graphic_$round\_$1.png\" USEMAP=\"#mapname$cont_html\" border=1><BR><BR>\n\n";
	    $cont_map++;
	}
        ##########################END graphic


	print HTML "Total number of reads in the last contig(s): $quant_reads<BR>\n";

	$/ = ">";
	open (CONTIGS, "fasta_cap_final.fasta.cap.contigs");
	<CONTIGS>;
	my $contig;
	while (<CONTIGS>){
	    chomp;
	    for my $contig_number_aux (keys %$contig_number){
		
		if ($_ =~ m/Contig$contig_number_aux/){
		    $contig .= ">$_";
		}
	    }
	}
	close (CONTIGS);
	

	my %contig_number = %$contig_number;
	my @header_original = @$header_original;
	
	my $number_seed_header = @header_original+0;



	for my $contig_number_aux (keys %$contig_number){
	    
	    if ($number_seed_header > 1){

		my $seed_name_aux = $contig_number{$contig_number_aux};
		$seed_name_aux =~ s/_//;

		my $name = "$seed_name_aux - final extended sequence - Contig$contig_number_aux";
		$contig =~ s/>Contig$contig_number_aux/<BR>>$name/;
	    }
	    else{
		my $name = "$header_original[0] - final extended sequence - Contig$contig_number_aux";
		$contig =~ s/>Contig$contig_number_aux/<BR>>$name/;
	    }
	}

	my $contig_html = $contig;
	$contig_html =~ s/\n/<BR>\n/g;

	open (HTML_CONTIGS, ">>seed_contigs.html");
	print HTML_CONTIGS "FINAL CONTIG(S)<BR>\n";
	print HTML_CONTIGS "\n<FONT FACE=\"Courier\">$contig_html</FONT>\n</body>\n</html>\n";
	close(HTML_CONTIGS);

    }
    close (HTML);
    
    return;
}



sub gif_generation {

    my $file = shift;
    my $contigs_numbers = shift;
    my $round = shift;
    my $seed = shift;
    my $final = shift;
    my $seed_type = shift;
    my $rep_seed_position = shift;


    my $cont = 0;
    my %begin_read;
    my $name;
    my $y;
    my $white;
    my $black;
    my $red;
    my $blue;
    my $contig_number;
    my $contig_length;
    my $cont_reads  = 0;
    my $contro_new_contig;

    my $control_contig = 0;

    my @map;
    my @head;
    my $map;
    my $head;
    my $post_seed;
    my $end_old = 0;
    my @lines;
    my @seed_lines;
    my $contig_number_rep_seed;
    my $contig_number_rep_seed_old;
    my $cont1 = 0;
    my $graphic_contig_number;

    open (ACE_FILES, "$file");

    while (<ACE_FILES>){
	if ($_ =~ m/CO\s+Contig(\d+)\s+(\d+)/){
	    $cont++;
	
	    $contig_number_rep_seed  = $1;

	    $contig_number = $1;
	    $contig_length = $2;


	    if ($control_contig == 1){
		if ($cont > 1 ){
		    
		    my $cont_aux = 0;
		    foreach my $aux (@lines){
			
			my @aux = split (/ /, $aux);
			
			my $ini = $aux[0];
			my $end = $aux[1];
			my $name = $aux[2];
			
			if ($round > 1){
			    if ($name =~ m/contig_\d+\_seed_program/){
				push @seed_lines, "$ini $end $name" ;
				delete $lines[$cont_aux];
			    }	
			}
			$cont_aux++;
		    }

		    my @lines_sort = sort {$a <=> $b} @lines;
		    
		    @lines = ();
		    
		    my $end_old = 0;
		    my $end_max;
		    my $cont_seed = 0;
		    my $control_print_seed = 0;
		    foreach my $aux (@lines_sort){
			
			my @aux = split (/ /, $aux);
			
			my $ini = $aux[0];
			my $end = $aux[1];
			my $name = $aux[2];
			
			if (($end_old < $ini) && ($end_old > 0) && ($ini > $end_max)){
			    
			    my @aux = split (/ /,$seed_lines[$cont_seed]);
			    $cont_seed++;
			    my $ini = $aux[0];
			    my $end = $aux[1];
			    my $name = $aux[2];
			    my ($map_aux,$head_aux) = lines ($ini,$end,$name,$round);
			    
			    $map .= $map_aux;
			    $head .= $head_aux;
			    $control_print_seed++;
			}
			if ($end > $end_max){
			    $end_max = $end;
			}
			$end_old = $end;
			
			my ($map_aux,$head_aux) = lines ($ini,$end,$name,$round);
			
			$map .= $map_aux;
			$head .= $head_aux;
		    }
		    
		    my $seed_aux = @seed_lines+0;
		    if ($control_print_seed < $seed_aux){
			
			my @aux = split (/ /,$seed_lines[$control_print_seed]);
			$cont_seed++;
			my $ini = $aux[0];
			my $end = $aux[1];
			my $name = $aux[2];
			my ($map_aux,$head_aux) = lines ($ini,$end,$name,$round);
			$map .= $map_aux;
			$head .= $head_aux;
		    }
		    
		    
		    push @map, $map;
		    push @head, $head;
		    
		    $map = "";
		    $head = "";
		    
		    print FIGURE $im->png;
		    close(FIGURE);
		    @seed_lines = ();
		}
	    }

	    $contro_new_contig = 1;
	}
	if ($_ =~ m/AF\s+(.*)\s+(\w)\s+([-\d]+)/){

	    $begin_read{$1} = $3;
	    $cont_reads++;
	}
	if ($_ =~m/RD\s+(\S+)/){

	    $name = $1;
	    
	    if ($contro_new_contig == 1){
		
		for my $aux (@$contigs_numbers){
		    
		    if ($aux == $contig_number){

			$graphic_contig_number .= "$contig_number ##";
			$post_seed = graph($contig_length,$contig_number,$cont_reads,$cont);
			
			#if ($seed_type == 1){
			if (($round == 1) || ($final == 1)){
			    my %rep_seed_position_aux = %$rep_seed_position;
			    
			    for (keys %rep_seed_position_aux){
				$aux =~ m/(\d+)_/;
				
				if ($_ == $1){
				    
				    my @rep_seed_aux = split (/\#\#/, $rep_seed_position_aux{$_});
				    push @lines, @rep_seed_aux;
				}
				
			    }
			}

			$control_contig = 1;
			last;
		    }
		    else{
			$control_contig = 0;
		    }
		}


		$contro_new_contig = 0;
		$cont_reads=0;
	    }
	}

	if ($control_contig == 1){
	    if ($_ =~m/QA\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/){
		my $ini = $begin_read{$name}+$1;
		my $end = $begin_read{$name}+$2;
		push @lines, "$ini $end $name" 
	    }
	}
    }
    close (ACE_FILES);
######################################
###################close ACE FILE
######################################


    my $cont_aux = 0;
    foreach my $aux (@lines){
	
	my @aux = split (/ /, $aux);
	
	my $ini = $aux[0];
	my $end = $aux[1];
	my $name = $aux[2];

	if ($round > 1){
	    if ($name =~ m/contig_\d+\_seed_program/){
		push @seed_lines, "$ini $end $name" ;
		delete $lines[$cont_aux];
	    }	
	}
	$cont_aux++;
    }

    my @lines_sort = sort {$a <=> $b} @lines;

    my $end_old = 0;
    my $end_max;
    my $cont_seed = 0;
    my $control_print_seed = 0;

    foreach my $aux (@lines_sort){

	my @aux = split (/ /, $aux);

	my $ini = $aux[0];
	my $end = $aux[1];
	my $name = $aux[2];

	if (($end_old < $ini) && ($end_old > 0) && ($ini > $end_max)){
	    
	    my @aux = split (/ /,$seed_lines[$cont_seed]);
	    $cont_seed++;
	    my $ini = $aux[0];
	    my $end = $aux[1];
	    my $name = $aux[2];
	    my ($map_aux,$head_aux) = lines ($ini,$end,$name,$round);
	    $map .= $map_aux;
	    $head .= $head_aux;
	    $control_print_seed++;
	}

	if ($end > $end_max){
	    $end_max = $end;
	}

	$end_old = $end;

	my ($map_aux,$head_aux) = lines ($ini,$end,$name,$round);

	$map .= $map_aux;
	$head .= $head_aux;
    }
    
    my $seed_aux = @seed_lines+0;
    if ($control_print_seed < $seed_aux){
	
	my @aux = split (/ /,$seed_lines[$control_print_seed]);
	$cont_seed++;
	my $ini = $aux[0];
	my $end = $aux[1];
	my $name = $aux[2];
	my ($map_aux,$head_aux) = lines ($ini,$end,$name,$round);
	$map .= $map_aux;
	$head .= $head_aux;
    }

    push @map, $map;
    push @head, $head;

    print FIGURE $im->png;
    close(FIGURE);

    return (\@map,\@head,$graphic_contig_number);

    sub graph{
	
	my $width_size = shift;
	my $cont_number = shift;
	my $cont_reads = shift;
	my $cont = shift;

	my $alt_gif = ($cont_reads*5)+60;

	$y = 35;
	
	my $tamline = $width_size;
	$width_size =  $width_size/10;
	
	open (FIGURE, ">image_dir/graphical_representation$cont.png") or die "Can't open image_dir/graphical_representation$cont.png";
	$im = new GD::Image($width_size+40,$alt_gif);
	
	$white = $im->colorAllocate(255,255,255);
	$black = $im->colorAllocate(0,0,0);
	$red = $im->colorAllocate(255,0,0);
	$blue = $im->colorAllocate(0,0,255);
	
        #create name of fasta
	$im->string( gdSmallFont, 0,0 ,"Contig$cont_number", $black );
	
	#create reference
	$im->string( gdSmallFont, 10,15,0, $black );
	$im->string( gdSmallFont, $width_size,15,$tamline, $black );
	
	#create mark lines
	$im->line(10,30,10,25,$black);
	$im->line($width_size+10,30,$width_size+10,25,$black);
	
    
	#create fasta line
	$im->line(10,30,$width_size+10,30,$black);
	$im->line(10,31,$width_size+10,31,$black);
	$im->line(10,32,$width_size+10,32,$black);
	
	my $grid_line = int($width_size/50);

	for (my $aux = 0; $aux <= $grid_line; $aux++){
	    my $start_point = (50*$aux)+10;
	    my $start_point_label = 500*$aux;
	    $im->dashedLine($start_point,$alt_gif,$start_point,28,$black);
	    next if ($aux == 0);
	    next if ($start_point > ($width_size)-20);
	    
	    $im->string( gdSmallFont,$start_point-5,15,$start_point_label, $black );
	}
	
	$im->dashedLine($width_size+10,$alt_gif,$width_size+10,28,$black);

	return ($alt_gif-15);
    }
    
    sub lines {
    
	my $ini = shift;
	my $end = shift;
	my $name = shift;
	my $round = shift;

	my $color = $blue;
	my $map;
	my $head;
	my $seed_control = 0;



	if ($round == 1){
	    open (SEED_NAME,"$seed_original_file_name");
	    my @seeds_name;
	    
	    while (<SEED_NAME>){
		if ($_ =~ m/^>(.*)/){
		    push @seeds_name, $1;
		}
	    }
	    close (SEED_NAME);
	  
	    for my $aux (@seeds_name){

		if ($aux =~ m/$name/){
		    $color = $red;
		    $seed_control = 1;
		}
	    }

	    ### Avoid printing the final contig in the report.html file
	    if ($name =~ m/contig_\d+_seed_program/){
		next;
	    }

	    if ($name =~ m/representative_seed/){
		$color = $red;
		$seed_control = 1;
	    }
	}
	else{
	    if ($name =~ m/contig_(\d+)\_seed_program_(.*)/){
		$color = $red;
		$seed_control = 1;
		$name = "$2 - extended sequence - Contig $1";
	    }	
	}

	my $ini_aux = ($ini/10)+10;
	my $end_aux = ($end/10)+10;
	
	$y = $y+5;

	$im->line($ini_aux, $y,   $end_aux, $y,   $color);
	$im->line($ini_aux, $y+1, $end_aux, $y+1, $color);
	$im->line($ini_aux, $y+2, $end_aux, $y+2, $color);

	$number++;
	my $ymais3 = $y+3;
	$map = "<area shape=\"rect\" COORDS=\"$ini_aux,$y,$end_aux,$ymais3\"   href=\"JavaScript:Rep$number()\">\n";

	$ini = $ini-1;
	$end = $end-1;


	$head = "<script language=\"JavaScript\">
function Rep$number(){ 
  var repeatwin$number = window.open(\"\",
			     \"jan$number\",
			     \"height=200,width=300,scrollbars=YES\");
  repeatwin$number.document.write(\"<html>\\n<head>\\n<title>$name</title>\\n</head>\\n<body>\\n<b><B>Read- $name<BR><B><BR>Coordinates in the contig consensus:<BR>Start- $ini<BR>\\n<B>End- $end<BR>\\n\");
  repeatwin$number.document.write(\"<FORM>\\n\");
  repeatwin$number.document.write(\'<INPUT TYPE=\"BUTTON\" VALUE=\"Close\" onClick=\"window.close();\">\');
  repeatwin$number.document.write(\"</FORM>\\n</body>\\n</html>\\n\");
}
</script>
";
	return ($map,$head);
    }
}


### Finishing process 
### Create the last file to run CAP3
### Delete some temporary files
### Close the HTML report files

sub finishing {

    my $db = shift;
    my $cap_param = shift;
    my $header_original = shift;
    my $round = shift;
    $round--;

    print "####   Last Round   \####\n"; ### "\" only for Emacs indentation
    system "echo \"####   Last Round   \####\" >>genseed.log"; ### "\" only for Emacs indentation

    create_final_file($db);
    
    system "cat contig_$round.fasta >>fasta_cap_final.fasta" if ($put_last_contig_last_assembly =~ /yes/);

    #### Mask vector bases with Cross_match
    cross_match("fasta_cap_final.fasta",$vec_db,$cross_param,"1") if ($vec_db);

    system "cap3 fasta_cap_final.fasta $cap_param >cap3_temp_genomic 2>>genseed.log";

    my $contig_number;
    my $contig;
    my $contig_seed_name = "";
    my $rep_seed_position;
    my $last_contig_control;

    ($contig_number,$contig,$rep_seed_position,$last_contig_control) = parser_cap3($db,$header_original,"fasta_cap_final.fasta.cap.ace","1","",$contig_length,$round,"",$seed_orientation_control,$original_seed_type);


    if ($last_contig_control == 1){
    

	my %contig_number = %$contig_number;
	my @header_original = @$header_original;
	my $number_seed_header = @header_original+0;
	
	
	
	report("last","1",$contig_number,$header_original,$contig,$original_seed_type,$rep_seed_position);
	
	open (FINAL_CONTIG, "fasta_cap_final.fasta.cap.contigs");
	open (OUT_FINAL_CONTIG, ">final_contigs.fasta");
	my $header_print;
	
	$/ = ">";
	<FINAL_CONTIG>;
	while (<FINAL_CONTIG>){
	    chomp;
	    for my $contig_number_aux (keys %{$contig_number}){
		
		$_ =~ m/Contig(\d+)/;
		
		
		if ($contig_number_aux == $1){
		    
		    $_ =~ /(Contig$contig_number_aux)(.*)/s;
		    
		    
		    my $header = $1;
		    my $sequence_aux = $2;
		    
		    if ($number_seed_header > 1){
			my $aux_name = $contig_number{$contig_number_aux};
			$aux_name =~ s/_//;
			$header_print = "$aux_name - final extended sequence - $header";
		    }
		    else{
			$header_print = "$header_original[0] - final extended sequence - $header";
		    }
		    
		    
		    print OUT_FINAL_CONTIG ">$header_print\n";
		    
		    $sequence_aux =~ s/\n//g;
		    
		    while ($sequence_aux){
			my $sequence_aux = substr ($sequence_aux,0,60,"");
			print OUT_FINAL_CONTIG "$sequence_aux\n";	    
		    }
		}
	    }
	}
	close(FINAL_CONTIG);
	close (OUT_FINAL_CONTIG);
	
    }
    else{

	print "\nSeed sequence was not found in the final contig(s). Aborting process.\n\n";
	system "echo \"Seed sequence was not found in the final contig(s). Aborting process.\" >>genseed.log"; ### "\" only for Emacs indentation
    }

    #####move CAP3 files to output directory
    system "mv fasta_cap_final.fasta fasta_dir/final_reads.fasta";
    system "mv fasta_cap_final.fasta.cap.ace CAP3_dir/final_reads.fasta.cap.ace";
    system "mv fasta_cap_final.fasta.cap.contigs CAP3_dir/final_reads.fasta.cap.contigs";
    system "mv fasta_cap_final.fasta.cap.contigs.links CAP3_dir/final_reads.fasta.cap.contigs.links";
    system "mv fasta_cap_final.fasta.cap.contigs.qual CAP3_dir/final_reads.fasta.cap.contigs.qual";
    system "mv fasta_cap_final.fasta.cap.info CAP3_dir/final_reads.fasta.cap.info";
    system "mv fasta_cap_final.fasta.cap.singlets CAP3_dir/final_reads.fasta.cap.singlets";


    ##### Remove temporary files
    system "rm blast_output cap_input* cap3_temp_genomic fasta_cap.fasta contig_*.fasta";
    if (-e "accepted_reads.fasta"){
	system "rm accepted_reads.fasta";
    }
    if (-e "formatdb.log"){
	system "rm formatdb.log";
    }

    if (-e "seed_fasta1.fasta"){
	system "rm seed_fasta1.fasta";
    }

    if (-e "temp_pseq.fasta"){
	system "rm temp_pseq.fasta";
    }

    if (-e "accepted_reads.fasta.screen"){
	system "rm accepted_reads.fasta.screen";
    }
    if (-e "accepted_reads.fasta.log"){
	system "rm accepted_reads.fasta.log";
    }
    if (-e "fasta_cap_final.fasta.screen"){
	system "rm fasta_cap_final.fasta.screen";
    }
    if (-e "fasta_cap_final.fasta.log"){
	system "rm fasta_cap_final.fasta.log";
    }
    if (-e "cross_temp_genomic"){
	system "rm cross_temp_genomic";
    }

    if (-e "one_read.fasta"){
	system "rm one_read.fasta";
    }

    if ($out_files =~ /simple/){
	if (-e "CAP3_dir/cap_1.ace"){
	    system "rm CAP3_dir/cap_*.ace";
	    system "rm fasta_dir/consensus_*.fasta";
	}
    }
    
    system "rm $seed_original_file_name";
    
    chdir "../";

    exit;
}

sub create_final_file{

    my $blastdb = shift;

    for (keys %contigs){

	my $name = "\"$_\"";
	next if ($name =~ m/contig_\d+_seed_program/);

	system "fastacmd -d $blastdb -s $name >>fasta_cap_final.fasta 2>>genseed.log";
	
    }
}

### Seed_verify subroutine - verify type of seed (DNA or protein) and if the file is composed
  # by a single or multiple sequences in FASTA format
sub seed_verify {

    my $seed = shift;
    my $seed_type = 0;
    my $counter_fasta = 0;
    my $protein_sequence = ""; 

    open (DATA, "$seed");

    while (<DATA>){
	if ($_ =~ m/^>/){ 
	    $counter_fasta++;
	    next;
	}
	
	if ($_ =~ m/([lvipfsyqderkhwm]+)/i){
	    $seed_type = 1;
	    chomp;
	    $protein_sequence .= $_;
	}
    }
    return ($seed_type,$counter_fasta);
}
