from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-k", help="Desired kmer size for output queries.pep file")
parser.add_argument("-i", help="Input protists.pep file")
parser.add_argument("-o", help="Desired name for output file")
args = parser.parse_args()

kmer = int(args.k)

with open(args.o,'w') as out:
        for i in SeqIO.parse(args.i,'fasta'):
                for j in range(0,len(i.seq),20):
                        out.write(">"+str(i.id)+"_"+str(j)+"\n"+str(i.seq[j:j+kmer])+"\n")
