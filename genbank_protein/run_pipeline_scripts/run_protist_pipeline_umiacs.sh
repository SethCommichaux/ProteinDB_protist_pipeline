#!/bin/sh
#SBATCH --time=120:00:00
#SBATCH --mem=50GB
#SBATCH --cpus=12
#SBATCH --qos=large
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --job-name="protist"

# Load modules and software paths into environment
#
module load biopython trimmomatic

diamond="/fs/cbcb-scratch/scommich/ProtistDB_protein/diamond"
kaiju="/fs/cbcb-scratch/scommich/ProtistDB_protein/kaiju/bin/" 
createDB="/fs/cbcb-scratch/scommich/ProtistDB_protein/genbank_protein/createDB_scripts/"
data="/fs/cbcb-scratch/scommich/ProtistDB_protein/genbank_protein/data/"
run_pipeline="/fs/cbcb-scratch/scommich/ProtistDB_protein/genbank_protein/run_pipeline_scripts/"
kaijuDB="/fs/cbcb-scratch/scommich/ProtistDB_protein/genbank_protein/data/binningDB.fasta.kaiju.fmi"
queryDB="/fs/cbcb-scratch/scommich/ProtistDB_protein/genbank_protein/data/queryDB"
adapters="/cbcb/sw/RedHat-7-x86_64/common/local/trimomatic/0.36/share/trimmomatic-0.36-3/adapters/TruSeq3-SE.fa"


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
reads_fastq=${i%.fastq}.trimmed.fastq

# Output file
#
out=$output/$(basename ${i%.fastq})
mkdir $out


# Trim/Filter raw reads
#
trimmomatic SE -threads 12 $i $reads_fastq ILLUMINACLIP:$adapters MAXINFO:120:0.2


# Run kaiju to query fastq reads against protein sequence binning databse (binningDB.fasta)
#
$kaiju/kaiju -t $data/nodes.dmp -f $kaijuDB -i $reads_fastq -o $out/kaiju -z 12 -m 9
#$kaiju/addTaxonNames -t $data/nodes.dmp -n $data/names.dmp -i $out/kaiju -o $out/kaiju.taxa -u -p

# Extract reads that aligned to binning database
#
python $run_pipeline/extract_kaiju_reads.py -k $out/kaiju -s $reads_fastq -o $out/kaiju.fasta
rm $out/kaiju

# Align binned reads, with Diamond, to queryDB
#
time $diamond blastx --db $queryDB --query $out/kaiju.fasta --threads 12 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen --id 60 --subject-cover 30 --out $out/kaiju.fasta.diamond

# Process diamond output
#
python $run_pipeline/subset_diamond_best_bitscore.py -d $out/kaiju.fasta.diamond -o $out/kaiju.fasta.diamond.subset
python $run_pipeline/process_diamond_output.py -d $out/kaiju.fasta.diamond.subset -f $out/kaiju.fasta -t $data/genbank_map_ncbi_taxonomy.txt -p $data/proteinID_map_length.txt

# Create files for producing sankey diagrams
#
python $run_pipeline/sankey_preprocess.py -n $data/nodes.dmp -f $data/fullnamelineage.dmp -t $out/kaiju.fasta.diamond.subset.results -o $out/

done
