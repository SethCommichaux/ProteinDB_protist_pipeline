import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-k", help="Kaiju output file")
parser.add_argument("-u", help="UniRef100 fasta file")
parser.add_argument("-o", help="Desired output fasta file")
args = parser.parse_args()


kaiju = {i.strip().split('\t')[1].split(' ')[0]:0 for i in open(args.k)}

with open(args.o,'w') as out:
	for i in SeqIO.parse(args.u,'fasta'):
		if str(i.id) in kaiju:
			out.write(">"+str(i.description)+"\n"+str(i.seq)+"\n")
