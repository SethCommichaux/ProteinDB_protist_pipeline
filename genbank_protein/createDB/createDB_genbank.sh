#!/bin/sh
#SBATCH --time=200:00:00
#SBATCH --mem=50GB
#SBATCH --cpus=12
#SBATCH --qos=large
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --job-name="get_proteins"


# Load modules and software paths into environment
#
module load Python2/common/2.7.9
diamond="/fs/cbcb-scratch/scommich/ProtistDB_protein/diamond"
kaiju="/fs/cbcb-scratch/scommich/ProtistDB_protein/kaiju/bin/" 
createDB="/fs/cbcb-scratch/scommich/ProtistDB_protein/createDB_scripts/"
data="/fs/cbcb-scratch/scommich/ProtistDB_protein/data/"
run_pipeline="/fs/cbcb-scratch/scommich/ProtistDB_protein/run_pipeline_scripts/"


# Change to data directory
#
mkdir $data
cd $data


# Download, decompress and concatenate genbank protein sequence files 
#
wget ftp://ftp.ncbi.nih.gov/ncbi-asn1/protein_fasta/gb*.fsa_aa.gz
gunzip *gz
cat *.fsa_aa > queryDB.pep
rm *fsa_aa

# Download, decompress NCBI taxonomy files
#
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
tar xvf new_taxdump.tar.gz
rm new_taxdump.tar.gz rankedlineage.dmp taxidlineage.dmp type* citations.dmp delnodes.dmp division.dmp gencode.dmp host.dmp merged.dmp


# Make queryDB
#
$diamond makedb --in queryDB.pep --db queryDB --threads 12


# Create mapping file that matches genbank fasta sequence IDs to NCBI taxonomy
#
python $createDB/map_gbkID_2_taxonomy.py -g queryDB.pep -t fullnamelineage.dmp -o genbank_map_ncbi_taxonomy.txt
python $createDB/protein_length_map.py -f queryDB.pep -o proteinID_map_length.txt


# Extract protist protein sequences from genbank protein sequence file
#
time python $createDB/filter_genbank_for_protist_taxa.py -g queryDB.pep -l fullnamelineage.dmp -o binningDB.pep


# Create binning database
#
$kaiju/mkbwt -o binningDB.pep.kaiju -n 12 -l 100000 binningDB.pep
$kaiju/mkfmi binningDB.pep.kaiju
rm binningDB.pep.kaiju.bwt binningDB.pep.kaiju.sa


# Download, decompress and process protist BUSCO genes
#
#wget https://busco.ezlab.org/datasets/alveolata_stramenophiles_ensembl.tar.gz
#wget https://busco.ezlab.org/datasets/protists_ensembl.tar.gz
#tar xvf alveolata_stramenophiles_ensembl.tar.gz
#tar xvf protists_ensembl.tar.gz
#cat alveolata_stramenophiles_ensembl/ancestral* protists_ensembl/ancestral* > protist_ancestral_busco_proteins.fasta
#rm -r alveolata_stramenophiles_ensembl* protists_ensembl*

# Identify protist proteins homologous to protist BUSCO genes
#
#$diamond makedb --in protist_ancestral_busco_proteins.fasta --db protist_ancestral_busco_proteins --threads 12
#$diamond blastp --db protist_ancestral_busco_proteins --query binningDB.pep --threads 12 --outfmt 6 --id 50 --subject-cover 50 --out protist_busco_candidates.txt --max-target-seqs 1
