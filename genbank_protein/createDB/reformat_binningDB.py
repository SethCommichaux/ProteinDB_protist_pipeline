import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-d", help="binning database; fasta format")
parser.add_argument("-m", help="genbankID to NCBI taxonomy mapping file")
parser.add_argument("-o", help="desired name for reformatted binning database; fasta format")
args = parser.parse_args()

genbankID2taxaID = {}

for h,i in enumerate(open(args.m)):
	if h != 0:
	        try:
        	    	genbankID2taxaID[i.strip().split('\t')[0]] = i.strip().split('\t')[2]
	        except:
        	       	continue


with open(args.o,'w') as out:
	for i in SeqIO.parse(args.d,'fasta'):
		d = str(i.description)
		id = str(i.id)
		s = str(i.seq)
		if id in genbankID2taxaID:
			out.write(">"+genbankID2taxaID[id]+"\n"+s+"\n")
