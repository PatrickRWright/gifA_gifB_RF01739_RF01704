#!/usr/bin/env perl

use strict;
use warnings;

my $ctab = $ARGV[0];
my $genome_list = $ARGV[1];

# grep header
my $orgs = `grep 'ClusterID' $ctab`;
chomp $orgs;

my @split_orgs = split(/\t/, $orgs);

# remove non organism entries
shift @split_orgs;
shift @split_orgs;
shift @split_orgs;
pop @split_orgs;
pop @split_orgs;
pop @split_orgs;
pop @split_orgs;
pop @split_orgs;
pop @split_orgs;

foreach(@split_orgs) {
    # grep org specfic lines from genomelist.txt
    my $curr_line = `grep -P '\t$_\t' $genome_list`;
    chomp $curr_line;
    my @split_line = split(/\t/, $curr_line);
    my $ref_seq_ids = $split_line[11];
    #print "$_ ##" . $split_line[11] . "\n";
    my @split_refids = split(/,/, $ref_seq_ids);
    foreach(@split_refids) {
        # download fasta record
        system "~/gifA_gifB_RF01739_RF01704/scripts/get_fasta_from_refid.pl -acc $_ -g $_.fasta";
        print "~/gifA_gifB_RF01739_RF01704/scripts/get_fasta_from_refid.pl -acc $_ -g $_.fasta\n";
    }
}


