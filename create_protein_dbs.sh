#!/bin/sh
#$ -N get_proteins
#$ -j y
#$ -pe mpi 12
#$ -cwd

# Load necessary modules
module load biopython 
module load kaiju 
module load diamond

# Create temporary file for downloading protein sequences from RefSeq
mkdir protein_seqs
cd protein_seqs

# Download all Eukaryote RefSeq protein sequences from NCBI ftp site
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/fungi/fungi.*.protein.faa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/invertebrate/invertebrate.*.protein.faa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.*.protein.faa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plant/plant.*.protein.faa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.*.protein.faa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/protozoa/protozoa.*.protein.faa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian.*.protein.faa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_other/vertebrate_other.*.protein.faa.gz

# Decompress all files
gunzip *gz

# Concatenate files
cat *faa > ../RefSeq_Eukaryote_proteins.pep

# Remove temporary file with downloaded sequence data
cd ..
rm -r protein_seqs/

# Download and decompress SwissProt bacterial protein sequences
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_bacteria.dat.gz
gunzip uniprot_sprot_bacteria.dat.gz

# Convert SwissProt protein sequence files to fasta format using biopython
python /lustre/projects/SethCommichaux/MTDNA/Protist_pipeline_scripts/swiss2fasta.py uniprot_sprot_bacteria.dat
rm uniprot_sprot_bacteria.dat

# This step is to filter out eukaryote sequences that are similar to bacterial proteins
# That might have been contamination or from horizontal gene transfer. Removal is designed
# To reduce the number of false positive hits that come from bacterial sequences mapping to
# Eukaryotic proteins.

# Create Diamond index from SwissProt bacterial sequences
# Align eukaryote proteins to bacterial proteins saving unaligned sequences
diamond makedb --in uniprot_sprot_bacteria.fasta --db uniprot_sprot_bacteria
diamond blastp --query RefSeq_Eukaryote_proteins.pep --db uniprot_sprot_bacteria --threads 12 --outfmt 6 --id 95 --query-cover 95 --sensitive --max-target-seqs 1 --un RefSeq_Eukaryote_proteins.pep.unaligned --out eukaryote2bacteria.txt

# Reformat fasta headers to be NCBI taxonomic identifiers
python get_ids.py -i RefSeq_Eukaryote_proteins.pep.unaligned -o RefSeq_Eukaryote_proteins.taxa.pep -taxa_lineage /lustre/projects/SethCommichaux/MTDNA/Protist_pipeline_DB/fullnamelineage.dmp

# Create diamond index of eukaryote refseq protein database
diamond makedb --in RefSeq_Eukaryote_proteins.taxa.pep --db RefSeq_Eukaryote_proteins.taxa

# Create kaiju index of eukaryote refseq protein database
mkbwt -o RefSeq_Eukaryote_proteins.taxa.pep.kaiju -n 12 -l 100000 RefSeq_Eukaryote_proteins.taxa.pep
mkfmi RefSeq_Eukaryote_proteins.taxa.pep.kaiju

# Delete intermediate kaiju files
rm RefSeq_Eukaryote_proteins.taxa.pep.kaiju.bwt RefSeq_Eukaryote_proteins.taxa.pep.kaiju.sa
