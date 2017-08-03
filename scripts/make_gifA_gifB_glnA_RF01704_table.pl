#!/usr/bin/env perl

use strict;
use warnings;

# builds the comparison table for gifA, gifB, glnA and RF01704

# grep header
my $orgs = `grep 'ClusterID' ~/gifA_gifB_RF01739_RF01704/cluster.tab`;
chomp $orgs;

# returns the kegg ids
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

my %refid_to_kegg = (); # NC_000911 -> syn
my %kegg_to_orgname = (); # syn -> Synechocystis sp. PCC 6803

foreach(@split_orgs) { # @split_orgs are kegg ids
    # grep org specfic lines from genomelist.txt
    my $curr_line = `grep -P '\t$_\t' ~/gifA_gifB_RF01739_RF01704/genomelist.txt`;
    chomp $curr_line;
    my @split_line = split(/\t/, $curr_line);
    my $organism = $split_line[0];
    $organism =~ s/;/_/g;
    my $kegg_id = $split_line[1];
    $kegg_to_orgname{$kegg_id} = $organism;
    my $ref_seq_ids = $split_line[11];
    my @split_refids = split(/,/, $ref_seq_ids);
    # there may be several RefSeq IDs
    foreach(@split_refids) {
        $refid_to_kegg{$_} = $kegg_id; 
    }
}

# test
#print $refid_to_kegg{'NC_000911'} . "\n";
#print $kegg_to_orgname{$refid_to_kegg{'NC_000911'}} . "\n";

# if an entry is present then the specific organism has a gifB homolog
my %gifB_hash = (); # kegg_id -> count

my $gifB_ctab = `cat ~/gifA_gifB_RF01739_RF01704/SLL1515_cluster.tab`;
chomp $gifB_ctab;
my @split_gifB = split(/\t+/, $gifB_ctab);
foreach(@split_gifB) {
    # this expression returns the kegg id
    if ($_ =~ m/(\w+):/) {
        # quantify the amount of matches
        my $match_count = `grep -o "$1:" ~/gifA_gifB_RF01739_RF01704/SLL1515_cluster.tab | wc -l`;
        chomp $match_count;
        $gifB_hash{$1} = $match_count;
    }
}

# if an entry is present then the specific organism has a gifA homolog
my %gifA_hash = (); # kegg_id -> count

my $gifA_ctab = `cat ~/gifA_gifB_RF01739_RF01704/SSL1911_cluster.tab`;
chomp $gifA_ctab;
my @split_gifA = split(/\t+/, $gifA_ctab);
foreach(@split_gifA) {
    # this expression returns the kegg id
    if ($_ =~ m/(\w+):/) {
        # quantify the amount of matches
        my $match_count = `grep -o "$1:" ~/gifA_gifB_RF01739_RF01704/SSL1911_cluster.tab | wc -l`;
        chomp $match_count;
        $gifA_hash{$1} = $match_count;
    }
}

# if an entry is present then the specific organism has a glnA homolog
my %glnA_hash = (); # kegg_id -> count
my $glnA_refids = `grep -oP "NC_\\d{6}|NZ_.{8}" ~/gifA_gifB_RF01739_RF01704/cmsearch_RF01739_HP.txt | tr '\n' ','`;
chomp $glnA_refids;
my @split_glnA_refids = split(/,/, $glnA_refids);
foreach(@split_glnA_refids) {
    $glnA_hash{$refid_to_kegg{$_}}++;
}

# RF01704
my %RF01704_hash = (); # kegg_id -> count
my $RF01704_refids = `grep -oP "NC_\\d{6}|NZ_.{8}" ~/gifA_gifB_RF01739_RF01704/cmsearch_RF01704_HP.txt | tr '\n' ','`;
chomp $RF01704_refids;
my @split_RF01704_refids = split(/,/, $RF01704_refids);
foreach(@split_RF01704_refids) {
    $RF01704_hash{$refid_to_kegg{$_}}++;
}

my $RF01739_pfamscan = $ARGV[0];

# glnA downstream peptides
my %glnA_pep_hash = (); # kegg_id -> list of domains (DUF4278,DUF772) 
open(MYDATA, $RF01739_pfamscan) or die("\nError: cannot open $RF01739_pfamscan at make_gifA_gifB_glnA_RF01704_table.pl\n\n");
    my @glnA_pfamscan_lines = <MYDATA>;
close MYDATA;

foreach (@glnA_pfamscan_lines) {

    chomp $_;
    # skip header
    next if ($_ =~ m/#/);

    # check for refseq ID
    if ($_ =~ m/(NC_\d{6}|NZ_.{8})/) {
        my $RefSeq_ID = $1;
        # print $RefSeq_ID . "\n";
        my @split_pfam_line = split(/\s+/, $_);
        my $domain = $split_pfam_line[6];
        if (exists $glnA_pep_hash{$refid_to_kegg{$RefSeq_ID}}) { # if it already exists the list is extended
            $glnA_pep_hash{$refid_to_kegg{$RefSeq_ID}} = $glnA_pep_hash{$refid_to_kegg{$RefSeq_ID}} . ",$domain";
        } else {
            $glnA_pep_hash{$refid_to_kegg{$RefSeq_ID}} = $domain;
        }
    }
}

my $RF01704_pfamscan = $ARGV[1];

# RF01704 downstream peptides
my %RF01704_pep_hash = (); # kegg_id -> list of domains (DUF4278,DUF772)
open(MYDATA, $RF01704_pfamscan) or die("\nError: cannot open $RF01704_pfamscan at make_gifA_gifB_glnA_RF01704_table.pl\n\n");
    my @RF01704_pfamscan_lines = <MYDATA>;
close MYDATA;

foreach (@RF01704_pfamscan_lines) {

    chomp $_;
    # skip header
    next if ($_ =~ m/#/);

    # check for refseq ID
    if ($_ =~ m/(NC_\d{6}|NZ_.{8})/) {
        my $RefSeq_ID = $1;
        my @split_pfam_line = split(/\s+/, $_);
        my $domain = $split_pfam_line[6];
        if (exists $RF01704_pep_hash{$refid_to_kegg{$RefSeq_ID}}) { # if it already exists the list is extended
            $RF01704_pep_hash{$refid_to_kegg{$RefSeq_ID}} = $RF01704_pep_hash{$refid_to_kegg{$RefSeq_ID}} . ",$domain";
        } else {
            $RF01704_pep_hash{$refid_to_kegg{$RefSeq_ID}} = $domain;
        }
    }
}

# print table header 
print "org_name;org_abbrev;replicon_RefSeq_ids;gifA;gifB;glnA_RF01739;glnA_RF01739_positions;glnA_RF01739_downstream_peptide;RF01704;RF01704_positions;RF01704_downstream_peptide\n";

foreach(@split_orgs) {

    my $print_line = "";
    # add org_name
    $print_line = $kegg_to_orgname{$_};
    # add kegg_id (abbrev)
    $print_line = $print_line . ";$_";

    # add RefSeq ids and prepare RNA family hit positions
    my $RefSeq_ID_list = "";
    my $RF01704_positions = "";
    my $RF01739_positions = "";
    foreach my $key (sort keys %refid_to_kegg) {
        # if the RefSeq id matches the current organism (abbrev)
        if ($refid_to_kegg{$key} eq $_) {
            $RefSeq_ID_list = $RefSeq_ID_list . "$key,";
            # prepare RNA family position output
            # RF01704
            my $RF01704_grep = `grep "$key" ~/gifA_gifB_RF01739_RF01704/cmsearch_RF01704_HP.txt | awk '{ print \$7 " " \$8 " " \$9 }' | tr '\n' ','`;
            chomp $RF01704_grep if ($RF01704_grep);
            $RF01704_positions = $RF01704_positions . "$key: " . $RF01704_grep if ($RF01704_grep);
            # RF01739
            my $RF01739_grep = `grep "$key" ~/gifA_gifB_RF01739_RF01704/cmsearch_RF01739_HP.txt | awk '{ print \$7 " " \$8 " " \$9 }' | tr '\n' ','`;
            chomp $RF01739_grep if ($RF01739_grep);
            $RF01739_positions = $RF01739_positions . "$key: " . $RF01739_grep if ($RF01739_grep);
        } 
    }

    # remove trailing comma
    chop $RF01704_positions;
    chop $RF01739_positions;
    chop $RefSeq_ID_list;
    
    $print_line = $print_line . ";$RefSeq_ID_list";
    # add gifA
    if (exists $gifA_hash{$_}) {
        $print_line = $print_line . ";$gifA_hash{$_}";
    } else {
        $print_line = $print_line . ";0";
    }
    # add gifB
    if (exists $gifB_hash{$_}) {
        $print_line = $print_line . ";$gifB_hash{$_}";
    } else {
        $print_line = $print_line . ";0";
    }
    # add glnA
    if (exists $glnA_hash{$_}) {
        $print_line = $print_line . ";$glnA_hash{$_}";
    } else {
        $print_line = $print_line . ";0";
    }
    # add glnA positions
    if ($RF01739_positions) {
        $print_line = $print_line . ";$RF01739_positions";
    } else {
        $print_line = $print_line . ";0";
    }
    # add glnA downstream peptide
    if (exists $glnA_pep_hash{$_}) {
        $print_line = $print_line . ";$glnA_pep_hash{$_}";
    } else {
        $print_line = $print_line . ";0";
    }
    # add RF01704
    if (exists $RF01704_hash{$_}) {
        $print_line = $print_line . ";$RF01704_hash{$_}";
    } else {
        $print_line = $print_line . ";0";
    }
    # add RF01704 positions
    if ($RF01704_positions) {
        $print_line = $print_line . ";$RF01704_positions";
    } else {
        $print_line = $print_line . ";0";
    }
    # add RF01704 downstream peptide
    if (exists $RF01704_pep_hash{$_}) {
        $print_line = $print_line . ";$RF01704_pep_hash{$_}";
    } else {
        $print_line = $print_line . ";0";
    }

    print $print_line . "\n";
}


