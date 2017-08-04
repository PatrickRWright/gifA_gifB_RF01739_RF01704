# directory ~/gifA_gifB_RF01739_RF01704
# these commands assume that you have cloned 
# gifA_gifB_RF01739_RF01704 into your home directory

# PfamScan can be downloaded from ftp://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/
# and can be set up according to the README given within the package (PfamScan.tar.gz v.1.6)
# for this set of commands it is assumed that there is a directory ~/gifA_gifB_RF01739_RF01704/PfamScan

# get cluster.tab from http://mbgd.genome.ad.jp/ 
# (params -o1 -HO -S -c60 -p0.5 -V0.6 -C80 -ao0.8 -ai0.95 -ne1 -EVAL=0.00001 -SCORE=60 Last update 2016/05/19)
# orgs: in cluster.tab header 
# // the file is already deposited in gifA_gifB_RF01739_RF01704

# download genomelist.txt from http://mbgd.genome.ad.jp/htbin/genomelist
# // the file is already deposited in gifA_gifB_RF01739_RF01704

# prepare fasta db
scripts/grep_orgs_from_ctab.pl cluster.tab genomelist.txt
mkdir fasta
mv *fasta fasta
# // a zipped copy of the fasta directory is deposited in gifA_gifB_RF01739_RF01704
#    thus you can alternatively choose to "unzip fasta.zip" 

# make one fasta db
cat fasta/*fasta > 60_cyanos.fasta
# remove spaces in header
sed -i 's/ /_/g' 60_cyanos.fasta
# remove comma
sed -i 's/,//g' 60_cyanos.fasta

# prepare cm for cmseach
mkdir Rfam
cd Rfam
# RF01739
wget http://rfam.xfam.org/family/RF01739/cm
mv cm RF01739_from_HP.cm

# RF01704
wget http://rfam.xfam.org/family/RF01704/cm
mv cm RF01704_from_HP.cm

# calibrate model
cmcalibrate RF01739_from_HP.cm
cmcalibrate RF01704_from_HP.cm
# // the calibrated files are deposited in Rfam directory 
cd ..

# search calibrated cm against fasta db (evalue cutoff 0.0001)
cmsearch --noali -E 0.0001 Rfam/RF01739_from_HP.cm 60_cyanos.fasta > cmsearch_RF01739_HP.txt
cmsearch --noali -E 0.0001 Rfam/RF01704_from_HP.cm 60_cyanos.fasta > cmsearch_RF01704_HP.txt
# // the cmsearch output files are deposited in gifA_gifB_RF01739_RF01704

# prepare overlap
grep "SSL1911" cluster.tab >  SSL1911_cluster.tab
grep "SLL1515" cluster.tab > SLL1515_cluster.tab

# make subsequence fasta
# transeq from EMBOSS:6.5.7.0
# RF01739
scripts/parse_region_for_pfamscan.pl 600 cmsearch_RF01739_HP.txt 
mv fwd.fa RF01739_fwd.fa
transeq -sequence RF01739_fwd.fa -outseq RF01739_fwd_translated.fa -clean -frame F
mv rev.fa RF01739_rev.fa
# -frame is correctly set to F here because the reverse complement is built in parse_region_for_pfamscan.pl
transeq -sequence RF01739_rev.fa -outseq RF01739_rev_translated.fa -clean -frame F
cat RF01739_fwd_translated.fa RF01739_rev_translated.fa > RF01739_all_translated.fa

# RF01704
scripts/parse_region_for_pfamscan.pl 600 cmsearch_RF01704_HP.txt
mv fwd.fa RF01704_fwd.fa
transeq -sequence RF01704_fwd.fa -outseq  RF01704_fwd_translated.fa -clean -frame F
mv rev.fa RF01704_rev.fa
# -frame is correctly set to F here because the reverse complement is built in parse_region_for_pfamscan.pl
transeq -sequence RF01704_rev.fa -outseq RF01704_rev_translated.fa -clean -frame F
cat RF01704_fwd_translated.fa RF01704_rev_translated.fa > RF01704_all_translated.fa

# pfam_scan
cd PfamScan
# in ~/gifA_gifB_glnA_RF01704_analysis/PfamScan
./pfam_scan.pl -fasta ~/gifA_gifB_RF01739_RF01704/RF01739_all_translated.fa -dir . > ~/gifA_gifB_RF01739_RF01704/RF01739_subseqs_pfam_scan.txt
./pfam_scan.pl -fasta ~/gifA_gifB_RF01739_RF01704/RF01704_all_translated.fa -dir . > ~/gifA_gifB_RF01739_RF01704/RF01704_subseqs_pfam_scan.txt
# // the pfam_scan.pl output files are deposited in gifA_gifB_RF01739_RF01704 
cd ..

# build overlap table
scripts/make_gifA_gifB_glnA_RF01704_table.pl RF01739_subseqs_pfam_scan.txt RF01704_subseqs_pfam_scan.txt > gifA_gifB_glnA_RF01704_overlap_and_downstream_peptides.csv
# // the gifA_gifB_glnA_RF01704_overlap_and_downstream_peptides.csv coocurence table is deposited in gifA_gifB_RF01739_RF01704

# get 16S sequences
mkdir 16s_retrieval
# create list of RefSeq IDs
grep -oP "NC_\d{6}|NZ_.{8}" 60_cyanos.fasta > 16s_retrieval/RefSeq_ID_list.txt
# get gbk records // in 16s_retrieval
cd 16s_retrieval
# this command will download many genbank files
while read line ; do ~/gifA_gifB_RF01739_RF01704/scripts/get_refseq_from_refid.pl -acc $line -g $line.gb ; done <RefSeq_ID_list.txt
# print 16S sequences
while read line ; do ~/gifA_gifB_RF01739_RF01704/scripts/parse_16s_from_gbk.pl $line ; done <RefSeq_ID_list.txt > 16S_sequences.fa
# manually remove redundant Microcoleus_sp._PCC_7113_plasmid_pMIC7113.02_complete_sequence._NC_019760
# this organism seems to have a 16S rRNA encoded on its NC_019760 plasmid
sed -i '/NC_019760/,+1 d' 16S_sequences.fa
# // a copy of 16S_sequences.fa is deposited in 16s_retrieval
cd ..

