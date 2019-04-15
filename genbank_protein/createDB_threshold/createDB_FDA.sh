#! /bin/sh
#$ -N DB_create
#$ -j y
#$ -pe mpi 12
#$ -cwd

# Load modules and software paths into environment
#
module load biopython
module load pandas

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
grep -i 'scientific name' names.dmp | awk -F "\t|\t" '{print $1 "\t" $3}' > names
awk -F "\t|\t" '{print $1 "\t" $5}' nodes.dmp > nodes
rm nodes.dmp names.dmp new_taxdump.tar.gz rankedlineage.dmp taxidlineage.dmp type* citations.dmp delnodes.dmp division.dmp gencode.dmp host.dmp merged.dmp


# Download uniprot idmapping file and UniRef100 protein sequence file
#
#wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
#wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
#gunzip *gz


# Extract protist protein sequences
#
python $createDB/splitUniRef100.py -g uniref100.fasta -t fullnamelineage.dmp


# Download protist and eukaryota Busco proteins
#
wget https://busco.ezlab.org/datasets/protists_ensembl.tar.gz
wget https://busco.ezlab.org/datasets/eukaryota_odb9.tar.gz
tar xvfz protists_ensembl.tar.gz
tar xvfz eukaryota_odb9.tar.gz
mv protists_ensembl/ancestral_variants ancestral
cat eukaryota_odb9/ancestral_variants >> ancestral
$diamond makedb --in ancestral --db ancestral --threads 12
$diamond blastp --query protists.pep --db ancestral --threads 12 --id 30 --query-cover 50 --outfmt 6 --out busco_diamond.txt
python $createDB/extract_busco_diamond.py busco_diamond.txt protists.pep
awk -F "Tax=" '{print $2}'  busco.fasta | awk -F " TaxID" '{print $1}' | sort | uniq -c | sort -n > busco_taxa_counts.txt
$diamond makedb --in busco.fasta --db busco --threads 12
rm -r protists_ensembl* eukaryota_odb9*


# create negative database
#
python $createDB/get_negative.py busco_diamond.txt uniref100.fasta


# Create annotation file for protists.pep
#
python $createDB/processDB.py -q protists.pep -u idmapping.dat -t fullnamelineage.dmp
#rm idmapping.dat


python $createDB/kmer_split.py -k 30 -s 5 -i negatives.pep -o negatives30_5.pep
python $createDB/kmer_split.py -k 60 -s 10 -i negatives.pep -o negatives60_10.pep
python $createDB/kmer_split.py -k 90 -s 15 -i negatives.pep -o negatives90_15.pep
python $createDB/kmer_split.py -k 120 -s 20 -i negatives.pep -o negatives120_20.pep
python $createDB/kmer_split.py -k 30 -s 5 -i busco.pep -o busco30_5.pep
python $createDB/kmer_split.py -k 60 -s 10 -i busco.pep -o busco60_10.pep
python $createDB/kmer_split.py -k 90 -s 15 -i busco.pep -o busco90_15.pep
python $createDB/kmer_split.py -k 120 -s 20 -i busco.pep -o busco120_20.pep


for i in negatives*_*.pep;
do
$diamond blastp --query $i --db busco --threads 12 --id 50 --query-cover 100 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen --out $i.diamond
done

for i in busco*_*.pep;
do
$diamond blastp --query $i --db busco --threads 12 --id 50 --query-cover 100 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen --out $i.diamond
done


# Find thresholds for protist proteins
#
python find_thresholds_metaphyler.py -nodes nodes -names names -f protist_functions.txt
