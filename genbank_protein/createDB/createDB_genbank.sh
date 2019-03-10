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
cat *.fsa_aa > genbank.pep
rm *fsa_aa


# Download, decompress NCBI taxonomy files
#
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
tar xvf new_taxdump.tar.gz
rm new_taxdump.tar.gz rankedlineage.dmp taxidlineage.dmp type* citations.dmp delnodes.dmp division.dmp gencode.dmp host.dmp merged.dmp


# Create mapping file that matches genbank fasta sequence IDs to NCBI taxonomy
#
python $createDB/processDB.py -g genbank.pep -t fullnamelineage.dmp
rm genbank.pep


# Make queryDB
#
$diamond makedb --in queryDB.pep --db queryDB --threads 12


# Create binning database
#
$kaiju/mkbwt -o binningDB.pep.kaiju -n 12 -l 100000 binningDB.pep
$kaiju/mkfmi binningDB.pep.kaiju
rm binningDB.pep.kaiju.bwt binningDB.pep.kaiju.sa

