import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-f", help="Fasta file to extract ids and sequence lengths from")
parser.add_argument("-o", help="path name for desired output mapping file")
args = parser.parse_args()


with open(args.o,'w') as out:
	for i in SeqIO.parse(args.f,'fasta'):
		id = str(i.id)
		l = str(len(i.seq))
		out.write(id+"\t"+l+"\n")
