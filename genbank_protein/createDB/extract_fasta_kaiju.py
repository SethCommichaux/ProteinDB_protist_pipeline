import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-k", help="kaiju output file; output format 6")
parser.add_argument("-s", help="kaiju query fasta file")
parser.add_argument("-o", help="path name for output fasta file of extracted sequences")
args = parser.parse_args()

extract_ids = {i.strip().split('\t')[1].split(' ')[0]:0 for i in open(args.k) if i[0] == "C"}

with open(args.o,'w') as out:
	for i in SeqIO.parse(args.s,'fasta'):
		d = str(i.description)
		id = str(i.id)
		s = str(i.seq)
		if id in extract_ids:
			out.write(">"+d+"\n"+s+"\n")
