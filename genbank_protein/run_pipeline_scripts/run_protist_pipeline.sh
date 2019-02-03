#!/bin/sh
#SBATCH --time=10:00:00
#SBATCH --mem=50GB
#SBATCH --cpus=12
#SBATCH --qos=large
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --job-name="protist"

# Load modules and software paths into environment
#
diamond="/fs/cbcb-scratch/scommich/ProtistDB_protein/diamond"
kaiju="/fs/cbcb-scratch/scommich/ProtistDB_protein/kaiju/bin/" 
NCBItaxonomy="/fs/cbcb-scratch/scommich/Protist/protist_db/"
kaijuDB="/fs/cbcb-scratch/scommich/ProtistDB_protein/genbank_protein/data/binningDB.fasta.kaiju.fmi"
diamondDB="/fs/cbcb-scratch/scommich/ProtistDB_protein/genbank_protein/data/queryDB"

# Input and output directories
# 
input="/fs/cbcb-scratch/scommich/ProtistDB_protein/simulated_data/random10RefSeqGenomes/"
output="/fs/cbcb-scratch/scommich/ProtistDB_protein/simulated_data/random10RefSeqGenomes/"

mkdir $output

################################################################
################################################################
################################################################

for i in $input/*fastq

do 

# Fastq file(s) to be analyzed
#
reads_fastq=$i

# Output file
#
out=$output/$(basename ${i%.fastq})
mkdir $out

# Run kaiju to query fastq reads against protein sequence binning databse (binningDB.fasta)
#
$kaiju/kaiju -t $NCBItaxonomy/nodes.dmp -f $kaijuDB -i $reads_fastq -o $out/kaiju -z 12
$kaiju/addTaxonNames -t $NCBItaxonomy/nodes.dmp -n $NCBItaxonomy/names.dmp -i $out/kaiju -o $out/kaiju.taxa -u -p

# Extract reads that aligned to binning database
#
python extract_kaiju_reads.py -k $out/kaiju -s $reads_fastq -o $out/kaiju.fasta
rm $out/kaiju

# Align binned reads, with Diamond, to queryDB
#
time $diamond blastx --db $diamondDB --query $out/kaiju.fasta --threads 12 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen --id 60 --subject-cover 30 --out $out/kaiju.fasta.diamond

# Process diamond output
#
python subset_diamond_best_bitscore.py -d $out/kaiju.fasta.diamond -o $out/kaiju.fasta.diamond.subset
python process_diamond_output.py -d $out/kaiju.fasta.diamond.subset -f $out/kaiju.fasta -t /fs/cbcb-scratch/scommich/ProtistDB_protein/genbank_protein/data/genbank_map_ncbi_taxonomy.txt -p /fs/cbcb-scratch/scommich/ProtistDB_protein/genbank_protein/data/proteinID_map_length.txt

# Create files for producing sankey diagrams
#
python sankey_preprocess.py -n $NCBItaxonomy/nodes.dmp -f $NCBItaxonomy/fullnamelineage.dmp -t $out/kaiju.fasta.diamond.subset.results -o $out/

done
