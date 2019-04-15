from Bio import SeqIO
import sys

busco = {i.strip().split('\t')[0]:0 for i in open(sys.argv[1])}

with open('busco.pep','w') as out:
        for i in SeqIO.parse(sys.argv[2],'fasta'):
                id = str(i.id)
                s = str(i.seq)
                d = str(i.description)
                if id in busco:
                        out.write(">"+d+"\n"+s+"\n")
