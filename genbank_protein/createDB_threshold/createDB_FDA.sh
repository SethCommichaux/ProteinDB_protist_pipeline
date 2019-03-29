#! /bin/sh
#$ -N DB_create
#$ -j y
#$ -pe mpi 12
#$ -cwd

# Load modules and software paths into environment
#
module load biopython

diamond="/nfs/sw/apps/diamond/diamond"
kaiju="/nfs/sw/apps/kaiju/"
createDB="/lustre/projects/SethCommichaux/ProteinDB_protist_pipeline/createDB_scripts/"
data="/lustre/projects/SethCommichaux/ProteinDB_protist_pipeline/data/"
run_pipeline="/lustre/projects/SethCommichaux/ProteinDB_protist_pipeline/run_pipeline_scripts/"


# Change to data directory
#
mkdir $data
cd $data


# Download, decompress NCBI taxonomy files
#
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
tar xvf new_taxdump.tar.gz
rm new_taxdump.tar.gz rankedlineage.dmp taxidlineage.dmp type* citations.dmp delnodes.dmp division.dmp gencode.dmp host.dmp merged.dmp


# Download uniprot idmapping file and UniRef100 protein sequence file
#
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
gunzip *gz


# Extract protist protein sequences
#
python $createDB/splitUniRef100.py -g uniref100.fasta -t fullnamelineage.dmp
$diamond makedb --in protists.pep --db protists --threads 12
$kaiju/mkbwt -o protists.pep.kaiju -n 12 -l 100000 protists.pep
$kaiju/mkfmi protists.pep.kaiju
rm protists.pep.kaiju.bwt protists.pep.kaiju.sa 


# Extract homologs of protist proteins
#
$kaiju/kaijup -f protists.pep -i uniref100.fasta -z 12 -m 9 | grep "^C" > protist_homologs
python $createDB/extract_kaiju.py -k protist_homologs -o protist_homologs.pep -u uniref100.fasta
$diamond makedb --in protist_homologs.pep --db protist_homologs --threads 12
# rm uniref100.fasta protist_homologs


# Get homology information for each protist protein
#
python $createDB/kmer_split.py -k 40 -i protists.pep -o queries.pep
$diamond blastp --query queries.pep --db protist_homologs --threads 12 --id 50 --query-cover 100 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen --out protist_homologs.txt
# rm queries.pep protist_homologs.pep protist_homologs.dmnd


# Create annotation file for protists.pep
#
python $createDB/processDB.py -q protists.pep -u idmapping.dat -t fullnamelineage.dmp
#rm idmapping.dat


# Find protein-specific thresholds
#
python $createDB/find_thresholds.py -d protist_homologs.txt -m protist_functions.txt

