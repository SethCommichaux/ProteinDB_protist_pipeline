#!/bin/sh
#$ -N proteinDB
#$ -j y
#$ -pe mpi 12
#$ -cwd


# Load modules and software paths into environment
#
module load java

diamond="/nfs/sw/apps/diamond/diamond"
kaiju="/nfs/sw/apps/kaiju/"
createDB="/lustre/projects/SethCommichaux/ProteinDB_protist_pipeline/createDB_scripts/"
protist_data="/lustre/projects/SethCommichaux/ProteinDB_protist_pipeline/data/"
run_pipeline="/lustre/projects/SethCommichaux/ProteinDB_protist_pipeline/run_pipeline_scripts/"
kaijuDB="/lustre/projects/SethCommichaux/ProteinDB_protist_pipeline/data/binningDB.fasta.kaiju.fmi"
queryDB="/lustre/projects/SethCommichaux/ProteinDB_protist_pipeline/data/queryDB"
trimmomatic="/lustre/projects/SethCommichaux/ProteinDB_protist_pipeline/data/Trimmomatic-0.38/trimmomatic-0.38.jar"
adapters="/lustre/projects/SethCommichaux/ProteinDB_protist_pipeline/data/Trimmomatic-0.38/adapters/TruSeq3-SE.fa"


# Input and output directories
#
input="/lustre/projects/SethCommichaux/Kratom/data/"
output="/lustre/projects/SethCommichaux/Kratom/data/"
mkdir $output

################################################################
################################################################
################################################################

for i in $input/*fastq

do

# Fastq file(s) to be analyzed
#
reads_fastq=$(basename ${i%.fastq}.trimmed.fastq)

# Output file
#
out=$output/$(basename ${i%.fastq})
mkdir $out


# Trim/Filter raw reads
#
java -jar $trimmomatic SE -threads 12 $i $out/$reads_fastq ILLUMINACLIP:$adapters:2:30:10 MAXINFO:120:0.2


# Run kaiju to query fastq reads against protein sequence binning databse (binningDB.fasta)
#
$kaiju/kaiju -t $protist_data/nodes.dmp -f $kaijuDB -i $out/$reads_fastq -o $out/kaiju -z 12 -m 9

# Extract reads that aligned to binning database
#
python $run_pipeline/extract_kaiju_reads.py -k $out/kaiju -s $out/$reads_fastq -o $out/kaiju.fasta
rm $out/kaiju

# Align binned reads, with Diamond, to queryDB
#
time $diamond blastx --evalue 0.00001 --db $queryDB --query $out/kaiju.fasta --threads 12 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen --out $out/kaiju.fasta.diamond

# Process diamond output
#
python $run_pipeline/subset_diamond_best_bitscore.py -d $out/kaiju.fasta.diamond -o $out/kaiju.fasta.diamond.subset
python $run_pipeline/process_diamond_output.py -num_reads 2 -num_prots 2 -d $out/kaiju.fasta.diamond.subset -f $out/kaiju.fasta -t $protist_data/genbank_map_ncbi_taxonomy.txt -p $protist_data/proteinID_map_length.txt

# Create files for producing sankey diagrams
#
python $run_pipeline/sankey_preprocess.py -nodes $protist_data/nodes.dmp -lineage $protist_data/fullnamelineage.dmp -results $out/kaiju.fasta.diamond.subset.results -o $out

# Summarize analysis
#
python $run_pipeline/summarize.py -rr $i -tr $out/$reads_fastq -kf $out/kaiju -d $out/kaiju.fasta.diamond -ds $out/kaiju.fasta.diamond.subset -t $out/kaiju.fasta.diamond.subset.results -o $out/summary_report.txt

# Clean up
#
rm $out/$reads_fastq $out/kaiju.fasta
# $out/kaiju.fasta.diamond $out/kaiju.fasta.diamond.subset

done
