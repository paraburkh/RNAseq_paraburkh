#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

# Check for command-line arguments
my ($infile, $ann_ids, $read_names, $chr_loc) = @ARGV;

# $infile: refers to the fasta file whose IDs want to be renamed
# $ann_ids: refers to a file containing the protein ID in the field 1 and UniProt annotation in the second field
# $read_names: has read names on the first column and protein ID in the second column (only include annotated reads)

die "Usage: $0 <fastA_file> <ann_ids> <read_names> <chr_location>\n" unless $infile && $ann_ids;
# The <ann_ids> should have a ">" character at the beginning
# The <read_names> should also have the protein ID

# Store protein IDs and annotation into a hash
# example input: ADO3_ARATH	protein:ADO3_ARATH gene_symbol:ADO3 description:Adagio protein 3
my %ids;
open my $ID_FILE, '<', $ann_ids or die "Can't open $ann_ids: $!";
while (my $line = <$ID_FILE>) {
    chomp $line;
    my @temp = split('\t', $line);
    (my $id  = $temp[0]) =~ s/>//;
    my $ann = $temp[1];
    $ids{$id} = $ann;
}
close $ID_FILE;


# Store read ids in a hash using the protein IDs as values
# example input: >NODE_10003_length_1072_cov_49.608609_g9068_i0.p1	UBC7_ARATH
my %read_ids;
open my $READ_IDS, "<", $read_names or die "Can't open $read_names: $!";
while (my $line = <$READ_IDS>) {
    chomp $line;
    my @temp = split('\t', $line);
    (my $read_id  = $temp[0]) =~ s/>//;
    my $prot_id = $temp[1];
    $read_ids{$read_id} = $prot_id;
}
close $READ_IDS;


# Add another hash including read names and chromosome location and starting coordinate
my %chr_loc;
open my $CHR_LOC, "<", $chr_loc or die "Can't open $chr_loc: $!";
while (my $line = <$CHR_LOC>) {
    chomp $line;
    my ($read_id, @rest) = split('\t', $line);

    $chr_loc{$read_id} = \@rest;
}
close $CHR_LOC;


# Process the FASTA file and filter sequences
my $seq_in = Bio::SeqIO->new(
    -format   => 'fasta',
    -file     => $infile,
    -alphabet => 'dna'
);

my $protein_id;
my $annotation;
my $counter = 0;
my $formatted_counter;
while (my $seq = $seq_in->next_seq) {
    $counter++;
    $formatted_counter = sprintf("%03d", $counter);
    my $target_id = $seq->id;

    my $pre_id = $chr_loc{$target_id}->[0] . "." .  $chr_loc{$target_id}->[1] . "." . "p" . $formatted_counter;
    print "$pre_id\n";

    if (exists $read_ids{$target_id}){
            $protein_id = $read_ids{$target_id};
            $annotation = $ids{$protein_id}


    } else {
            $annotation = 'Unknown protein'
    }
    my $new_id = ">" . $pre_id . " " . $target_id . " " . $annotation; 
    #print $new_id, "\n", $seq->seq, "\n";
