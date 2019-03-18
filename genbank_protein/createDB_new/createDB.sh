#!/bin/sh
#SBATCH --time=200:00:00
#SBATCH --mem=80GB
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


# Create binningDB and nonprotistDB
#
python $createDB/splitUniRef100.py -g uniref100.fasta -t fullnamelineage.dmp
$kaiju/mkbwt -o binningDB.pep.kaiju -n 12 -l 100000 binningDB.pep
$kaiju/mkfmi binningDB.pep.kaiju
rm binningDB.pep.kaiju.bwt binningDB.pep.kaiju.sa uniref100.fasta


# Look for nonprotist proteins with 9-mer match to protist protein; combine with binningDB to create queryDB
#
$kaiju/kaijup -f binningDB.pep.kaiju.fmi -i nonprotists.pep -z 12 -m 9 | grep "^C" > protist_homologs.txt
rm nonprotists.pep


# Make queryDB
#
cp binningDB.pep queryDB.pep
python $createDB/extractKaiju.py -k protist_homologs.txt -np nonprotists.pep
$diamond makedb --in queryDB.pep --db queryDB --threads 12


# Create annotation file for queryDB
#
python $createDB/processDB.py -q queryDB.pep -u idmapping.dat -t fullnamelineage.dmp
#rm UniRef100_mappings.txt
