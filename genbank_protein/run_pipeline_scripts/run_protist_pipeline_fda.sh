#!/bin/sh
#$ -N proteinDB
#$ -j y
#$ -pe mpi 12
#$ -cwd


# Load modules and software paths into environment
#
module load java
module load biopython

diamond="/nfs/sw/apps/diamond/diamond"
kaiju="/nfs/sw/apps/kaiju/"
createDB="/lustre/projects/SethCommichaux/ProteinDB_protist_pipeline/createDB_scripts/"
protist_data="/lustre/projects/SethCommichaux/ProteinDB_protist_pipeline/data/"
run_pipeline="/lustre/projects/SethCommichaux/ProteinDB_protist_pipeline/run_pipeline_scripts/"
kaijuDB="/lustre/projects/SethCommichaux/ProteinDB_protist_pipeline/data/binningDB.fasta.kaiju.fmi"
queryDB="/lustre/projects/SethCommichaux/ProteinDB_protist_pipeline/data/queryDB"
trimmomatic="/lustre/projects/SethCommichaux/ProteinDB_protist_pipeline/Trimmomatic-0.38/trimmomatic-0.38.jar"
adapters="/lustre/projects/SethCommichaux/ProteinDB_protist_pipeline/Trimmomatic-0.38/adapters/TruSeq3-SE.fa"


# Input and output directories
#
input=""
output=""
mkdir $output

################################################################
################################################################
################################################################

for i in $input/*fastq

do

# Fastq file(s) to be analyzed
#
read_count=`grep -c "^+" $i`
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
$kaiju/kaijux -f $kaijuDB -i $out/$reads_fastq -z 12 -m 9 | grep "^C" > $out/kaiju

# Extract reads that aligned to binning database
#
python $run_pipeline/extract_kaiju_reads.py -k $out/kaiju -s $out/$reads_fastq -o $out/kaiju.fasta

# Align binned reads, with Diamond, to queryDB
#
time $diamond blastx --evalue 0.00001 --db $queryDB --query $out/kaiju.fasta --threads 12 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen --out $out/kaiju.fasta.diamond

# Process diamond output
#
python $run_pipeline/subset_diamond_best_bitscore.py -d $out/kaiju.fasta.diamond -o $out/kaiju.fasta.diamond.subset
python $run_pipeline/process_diamond_output.py -ANI 60 -num_reads 2 -num_prots 2 -d $out/kaiju.fasta.diamond.subset -f $read_count -t $protist_data/queryDB_functions.txt

# Create files for producing sankey diagrams
#
python $run_pipeline/sankey_preprocess.py -nodes $protist_data/nodes.dmp -lineage $protist_data/fullnamelineage.dmp -results $out/kaiju.fasta.diamond.subset.results -o $out

# Summarize analysis
#
echo "raw_fastq: " > $out/summary_report.txt; echo $read_count >> $out/summary_report.txt;
echo "trimmed_fastq: " >> $out/summary_report.txt; grep -c "^+" $out/$reads_fastq >> $out/summary_report.txt; 
echo "kaiju_fasta: " >> $out/summary_report.txt; grep -c "^>" $out/kaiju.fasta >> $out/summary_report.txt; 
echo "diamond_aligned: " >> $out/summary_report.txt; awk -F "\t" '{print $1}' $out/kaiju.fasta.diamond | sort | uniq | wc -l >> $out/summary_report.txt; 
echo "subset_diamond_aligned: " >> $out/summary_report.txt; awk -F "\t" '{print $1}' $out/kaiju.fasta.diamond.subset | sort | uniq | wc -l >> $out/summary_report.txt; 

# Clean up
#
rm $out/$reads_fastq $out/kaiju.fasta $out/kaiju
# $out/kaiju.fasta.diamond $out/kaiju.fasta.diamond.subset

done

