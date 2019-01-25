#!/bin/sh
#SBATCH --time=10:00:00
#SBATCH --mem=50GB
#SBATCH --cpus=12
#SBATCH --qos=normal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --job-name="protist"


# Load modules and software paths into environment
#
diamond="/fs/cbcb-scratch/scommich/ProtistDB_protein/diamond"
kaiju="/fs/cbcb-scratch/scommich/ProtistDB_protein/kaiju/bin/" 
NCBItaxonomy="/fs/cbcb-scratch/scommich/Protist/protist_db/"
kaijuDB="/fs/cbcb-scratch/scommich/ProtistDB_protein/genbank_protein/data/binningDB.fasta.kaiju.fmi"
diamondDB="/fs/cbcb-scratch/scommich/ProtistDB_protein/genbank_protein/data/queryDB"

# Fastq file(s) to be analyzed
#
reads_fastq="/fs/cbcb-scratch/scommich/ProtistDB_protein/simulated_data/random10RefSeqGenomes/10genomes.fastq"

# Output directory
#
output=${reads_fastq%.fastq}
mkdir $output

# Run kaiju to query fastq reads against protein sequence binning databse (binningDB.fasta)
#
$kaiju/kaiju -t $NCBItaxonomy/nodes.dmp -f $kaijuDB -i $reads_fastq -o $output/kaiju -z 12
$kaiju/addTaxonNames -t $NCBItaxonomy/nodes.dmp -n $NCBItaxonomy/names.dmp -i $output/kaiju -o $output/kaiju.taxa -u -p

# Extract reads that aligned to binning database
#
python extract_kaiju_reads.py -k $output/kaiju -s $reads_fastq -o $output/kaiju.fasta

# Align binned reads, with Diamond, to queryDB
#
time $diamond blastx --db $diamondDB --query $output/kaiju.fasta --threads 12 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen --id 80 --subject-cover 50 --out $output/kaiju.fasta.diamond

# Process diamond output
#
python process_diamond_output.py -d $output/kaiju.fasta.diamond -f $output/kaiju.fasta -t /fs/cbcb-scratch/scommich/ProtistDB_protein/genbank_protein/data/genbank_map_ncbi_taxonomy.txt

# Create files for producing sankey diagrams
#
python sankey_preprocess.py $output/kaiju.fasta.diamond.results
