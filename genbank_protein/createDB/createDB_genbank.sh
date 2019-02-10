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
createDB="/fs/cbcb-scratch/scommich/ProtistDB_protein/genbank_protein/createDB_scripts/"
data="/fs/cbcb-scratch/scommich/ProtistDB_protein/genbank_protein/data/"
run_pipeline="/fs/cbcb-scratch/scommich/ProtistDB_protein/genbank_protein/run_pipeline_scripts/"


# Change to data directory
#
cd $data


# Download, decompress and concatenate genbank protein sequence files
#
wget ftp://ftp.ncbi.nih.gov/ncbi-asn1/protein_fasta/gb*.fsa_aa.gz
gunzip *gz
cat *.fsa_aa > genbank.pep
rm *fsa_aa


# Download, decompress NCBI taxonomy files
#
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
tar xvf new_taxdump.tar.gz
rm new_taxdump.tar.gz rankedlineage.dmp taxidlineage.dmp type* citations.dmp delnodes.dmp division.dmp gencode.dmp host.dmp merged.dmp


# Extract protist protein sequences from genbank protein sequence file
#
time python $createDB/filter_genbank_for_protist_taxa.py -g genbank.pep -l fullnamelineage.dmp -o genbank_protists.pep
rm genbank.pep


# Create Kaiju DB and look for non-protist protein homologs of protist proteins
#
$kaiju/mkbwt -o genbank_protists.pep.kaiju -n 12 -l 100000 genbank_protists.pep
$kaiju/mkfmi genbank_protists.pep.kaiju
rm genbank_protists.pep.kaiju.bwt genbank_protists.pep.kaiju.sa
$kaiju/kaijup -f genbank_protists.pep.kaiju.fmi -i genbank.pepwithout_protists -o eukaryote2bacteria.txt -z 12


# Add non-protist protein homologs to protist proteins. Then create querying database
#
python $createDB/extract_fasta_kaiju.py -k eukaryote2bacteria.txt -s genbank.pepwithout_protists -o nonprotists.pep
cat nonprotists.pep genbank_protists.pep > queryDB.fasta
$diamond makedb --in queryDB.fasta --db queryDB --threads 12


# Create mapping file that matches genbank fasta sequence IDs to NCBI taxonomy
#
python $createDB/map_gbkID_2_taxonomy.py -g queryDB.fasta -t fullnamelineage.dmp -o genbank_map_ncbi_taxonomy.txt
python $createDB/protein_length_map.py -f queryDB.fasta -o proteinID_map_length.txt


# Create binning database
#
python $createDB/reformat_binningDB.py -d genbank_protists.pep -m genbank_map_ncbi_taxonomy.txt -o binningDB.fasta
$kaiju/mkbwt -o binningDB.fasta.kaiju -n 12 -l 100000 binningDB.fasta
$kaiju/mkfmi binningDB.fasta.kaiju
rm binningDB.fasta.kaiju.bwt binningDB.fasta.kaiju.sa


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


# Clean-up data directory
#
rm -r genbank.pep alveolata_stramenophiles_ensembl* protists_ensembl* eukaryote2bacteria.txt genbank.pepwithout_protists genbank_protists.dmnd genbank_protists.pep
