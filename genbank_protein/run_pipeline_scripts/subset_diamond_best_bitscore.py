import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d", type=str,
                        help="Input diamond results file (--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen)")
parser.add_argument("-o", type=str,
                        help="Desired output path for subsetted Diamond output file")
args = parser.parse_args()

read2bitscore = {}

with open(args.o,'w') as out:
	for i in open(args.d):
		tmp = i.strip().split('\t')
		read,bitscore = tmp[0],float(tmp[11])
		if read not in read2bitscore: 
			read2bitscore[read] = bitscore
			out.write(i)
		else:
			if bitscore == read2bitscore[read]:
				out.write(i)
				
