#!/usr/bin/env perl

use strict;
use warnings;

use Bio::SeqIO;

# this script retrieves the downstream region of RNA motifs and
# translates them into AA sequences to pass on to pfam_scan.pl 

my $downstream_length = $ARGV[0];
my $input_table = $ARGV[1];

open(MYDATA, $input_table) or die("\nError: Can't open file $input_table at parse_region_for_pfamscan.pl\n\n");
    my @intable_lines = <MYDATA>;
close(MYDATA);

open WRITEFWD, ">", "fwd.fa";
open WRITEREV, ">", "rev.fa";

my $c = 1;

foreach (@intable_lines) {

    chomp $_;
    # only read table
    next unless ($_ =~ m/\!/);
    # splitline
    my @split_line = split(/\s+/, $_);
    my $full_name = $split_line[6];
    my $start = $split_line[7];
    my $stop = $split_line[8];
    my $strand = $split_line[9];
   
    my $ID = "";
    # get refseq id
    if ($_ =~ m/(NC_\d{6}|NZ_.{8})/) {
        $ID = $1;
    }
     
    #print $ID . " start:$start stop:$stop $strand\n";

    my $curr_fasta = "fasta/$ID.fasta";
    my $seqio = Bio::SeqIO->new(-file => $curr_fasta, -format => "fasta");

    my $curr_seq = "";

    while (my $seq = $seqio->next_seq) {
        $curr_seq = $seq->seq();
    }

    # print $ID .": " . length($curr_seq) . "\n";

    my $fasta_header = ">$full_name" . "_$c\n";

    # always anchored on $stop

    my $subseq = "";

    if ($strand eq '-') {
        $subseq = reverseComplement(substr($curr_seq, (($stop-1)-$downstream_length), $downstream_length));
   #     print length($subseq) . "\n";
        print WRITEREV $fasta_header . $subseq . "\n";
    }
 
    if ($strand eq '+') {
        $subseq = substr($curr_seq, ($stop+1), $downstream_length);
   #     print length($subseq) . "\n";
        print WRITEFWD $fasta_header . $subseq . "\n";
    }
    $c++;

}

# get reverse complement of a nucleotide sequence
# 
sub reverseComplement{
    my ($seq) = @_;
    my $revcomp = reverse($seq);
    $revcomp =~ tr/ACGTUacgtu/TGCAAtgcaa/;
    return $revcomp;
}


