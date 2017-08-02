#!/usr/bin/env perl

use strict;
use warnings;

use Bio::SeqIO;

# parses 16s rRNA sequences from genbank files

my $RefID = $ARGV[0];
chomp $RefID;

my $seqin = Bio::SeqIO->new( -format => 'genbank', -file => "$RefID.gb");

while( (my $seq = $seqin->next_seq()) ) {
    foreach my $sf ( $seq->get_SeqFeatures() ) {
        if( $sf->primary_tag eq 'rRNA' ) {
            my $product = "";
            if ($sf->has_tag("product")) {
                my @productlist = $sf->get_tag_values("product");
                $product = $productlist[0];
            }

            if ($product =~ m/16S/i or $product =~ m/Small subunit ribosomal RNA/i) {
                # get organism description 
                my $organism = $seq->desc;
                $organism =~ s/ /_/g;
                $organism =~ s/,//g;
                # print fasta
                print ">$organism" . "_$RefID\n";
                print $sf->spliced_seq->seq, "\n";
                last;
            }

        }
    }
}


