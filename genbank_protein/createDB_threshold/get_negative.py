import sys
import random
from Bio import SeqIO

busco = {i.strip().split('\t')[0] for i in open(sys.argv[1])} # busco diamond file
uniref100 = sys.argv[2]

with open('negatives.pep','w') as out:
        for i in SeqIO.parse(uniref100,'fasta'):
                id,s,d = str(i.id),str(i.seq),str(i.description)
                if id not in busco:
                        x = random.randint(0,5) # 20% of sequences output
                        if x == 0:
                                out.write(">"+d+"\n"+s+"\n")
