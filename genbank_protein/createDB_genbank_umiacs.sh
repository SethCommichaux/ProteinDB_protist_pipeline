#!/bin/sh
#SBATCH --time=200:00:00
#SBATCH --mem=120GB
#SBATCH --cpus=12
#SBATCH --qos=large
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --job-name="get_proteins"


# Load modules and software paths into environment
#
module load Python2/common/2.7.9
diamond="/fs/cbcb-scratch/scommich/ProtistDB_protein/diamond"
kaiju="/fs/cbcb-scratch/scommich/ProtistDB_protein/kaiju/bin/" 


# Download, decompress and concatenate genbank protein sequence files
#
wget ftp://ftp.ncbi.nih.gov/ncbi-asn1/protein_fasta/gb*.fsa_aa.gz
gunzip *gz
cat *.fsa_aa > genbank.pep


# Create mapping file that matches genbank fasta sequence IDs to NCBI taxonomy
#
python map_gbkID_2_taxonomy.py -g genbank.pep -t ../../Protist/protist_db/fullnamelineage.dmp -o genbank_map_ncbi_taxonomy.txt
 

# Extract protist protein sequences from genbank protein sequence file 
#
time python filter_genbank_for_protist_taxa.py -g genbank.pep -l /fs/cbcb-scratch/scommich/Protist/protist_db/fullnamelineage.dmp -o genbank_protists.pep
rm genbank.pep


# Look for non-protist protein homologs of protist proteins
#
$diamond makedb --in genbank_protists.pep --db genbank_protists --threads 12
$diamond blastp --db genbank_protists --query genbank.pepwithout_protists --threads 12 --outfmt 6 --id 90 --subject-cover 50 --max-target-seqs 1 --out eukaryote2bacteria.txt


# Create binning database
#
cp genbank_protists.pep binningDB.fasta
$kaiju/mkbwt -o binningDB.fasta.kaiju -n 12 -l 100000 binningDB.fasta
$kaiju/mkfmi binningDB.fasta.kaiju
rm binningDB.fasta.kaiju.bwt binningDB.fasta.kaiju.sa binningDB.fasta


# Add non-protist protein homologs to protist proteins. Then create querying database
#
python extract_fasta_diamond.py -d eukaryote2bacteria.txt -s genbank.pepwithout_protists -o nonprotists.pep
cat nonprotists.pep genbank_protists.pep > queryDB.fasta
$diamond makedb --in queryDB.fasta --db queryDB --threads 12


# Download, decompress and process protist BUSCO genes
#
wget https://busco.ezlab.org/datasets/alveolata_stramenophiles_ensembl.tar.gz
wget https://busco.ezlab.org/datasets/protists_ensembl.tar.gz
tar xvf alveolata_stramenophiles_ensembl.tar.gz
tar xvf protists_ensembl.tar.gz
cat alveolata_stramenophiles_ensembl/ancestral* protists_ensembl/ancestral* > protist_ancestral_busco_proteins.fasta


# Identify protist proteins homologous to protist BUSCO genes
#
$diamond makedb --in protist_ancestral_busco_proteins.fasta --db protist_ancestral_busco_proteins --threads 12
$diamond blastp --db protist_ancestral_busco_proteins --query genbank_protists.pep --threads 12 --outfmt 6 --id 50 --subject-cover 50 --out protist_busco_candidates.txt --max-target-seqs 1
