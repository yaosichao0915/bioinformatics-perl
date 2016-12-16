#!/usr/bin/perl -w

=head1 NAME

codonalign.pl - align a coding sequence using a protein alignment as reference

=head1 SYNOPSIS

codonalign.pl aln_dir *.cut

=head1 DESCRIPTION

codonalign.pl align each fasta file using the corresponding 
clustalw .aln file in aln_dir as reference.

In the output alignment in fasta format, each aminoacid have been 
replaced by the corresponding three nucleotides from the coding sequence 
and each gap is replaced by three gaps.

The output files have the same name as the input files but with 
the .dna_aln extension. They are put in the current directory.

=head1 CAVEATS

The coding sequence file and alignment file must contains exactly 
the same sequences in the same order. There is no verification that 
the sequences really correspond.

A file having the same name as an output file will be overwritten.

=head1 AUTHOR

Based on align_on_codons.pl from bioperl and 
modified heavily by Jules Gagnon <eonwe@users.sourceforge.net>

=cut

use strict;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use File::Basename;

@ARGV >= 2 || die "Usage: codonalign.pl aln_dir *.cut\n";
my $CODONSIZE = 3;
my $alnpath   = shift;

while ( my $f = shift ) {

    my ( $name, $ext ) = split( /\./, basename($f) );

    my $seqio = new Bio::SeqIO(
        '-format' => 'fasta',
        '-file'   => $f
    );
    my $alignin = Bio::AlignIO->new(
        '-file'   => "$alnpath/$name.aln",
        '-format' => 'fasta'
    );
    my $alignout = new Bio::AlignIO(
        '-format' => 'fasta',
        '-file'   => ">$name.codon_aln"
    );

    my $aln      = $alignin->next_aln;
    my $alnlen   = $aln->length;
    my $dnaalign = new Bio::SimpleAlign;
    my $seqorder = 0;
    my @nucseqs;

    while ( my $seq = $seqio->next_seq ) {
        push @nucseqs, $seq;
        $aln->get_seq_by_pos( scalar(@nucseqs) )->id eq $seq->id
          || die "Sequence order mismatch "
          . $aln->get_seq_by_pos( scalar(@nucseqs) )->id . " vs "
          . $seq->id . "\n";
    }

    foreach my $seq ( $aln->each_seq ) {

        my $newseq;
        my $ngap = 0;
        my $end;
        foreach my $pos ( 1 .. $alnlen ) {
            my $loc = $seq->subseq( $pos, $pos );
            my $dna;
            if ( $loc eq '-' ) { $dna = '---'; $ngap++; }
            else {
                my $start = ( ( $pos - 1 - $ngap ) * $CODONSIZE ) + 1;
                $end = $start + $CODONSIZE - 1;
                die "Sequence length doesn't match for "
                  . $seq->id . " in "
                  . $name . "\n"
                  if ( $start > $nucseqs[$seqorder]->length
                    || $end > $nucseqs[$seqorder]->length );
                $dna = $nucseqs[$seqorder]->subseq( $start, $end );
            }
            $newseq .= $dna;
        }
        $end == $nucseqs[$seqorder]->length
          || die "Sequence length doesn't match for "
          . $seq->id . " in "
          . $name . "\n";
        my $newdna = new Bio::LocatableSeq(
            -display_id => $seq->id,
            -start      => ( ( $seq->start - 1 ) * $CODONSIZE ) + 1,
            -end        => $seq->end * $CODONSIZE,
            -strand     => $seq->strand,
            -seq        => $newseq
        );

        $dnaalign->add_seq($newdna);
        $seqorder++;
    }
    $dnaalign->set_displayname_flat;
    $alignout->write_aln($dnaalign);
}
