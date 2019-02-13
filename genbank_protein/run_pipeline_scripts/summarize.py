import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-rr", type=str,
                        help="Raw reads file")
parser.add_argument("-tr", type=str,
                        help="Trimmed reads file")
parser.add_argument("-kf", type=str,
                        help="Kaiju fasta file")
parser.add_argument("-d", type=str,
                        help="Diamond file")
parser.add_argument("-ds", type=str,
                        help="Diamond subset file")
parser.add_argument("-t", type=str,
                        help="Taxonomic results file")
parser.add_argument("-o", type=str,
                        help="Desired name and path for output file")
args = parser.parse_args()

rr,tr,kf = 0,0,0

for i in open(args.rr):
        if i[0] == "+": rr +=1

for i in open(args.tr):
        if i[0] == "+": tr +=1

for i in open(args.kf):
        if i[0] == ">": kf += 1

d = len({i.strip().split('\t')[0] for i in open(args.d)})

ds = len({i.strip().split('\t')[0] for i in open(args.ds)})

for i in open(args.t):
        if i.strip().split('\t')[0] == "CELLULAR ORGANISMS;EUKARYOTA;":
                t = [i.strip().split('\t')[1],i.strip().split('\t')[3],i.strip().split('\t')[4]]

with open(args.o,'w') as out:
        out.write("Total raw reads: "+str(rr)+'\n')
        out.write("Read count after trimming: "+str(tr)+'\n')
        out.write("Read count after Kaiju: "+str(kf)+'\n')
        out.write("Read count after diamond: "+str(d)+'\n')
        out.write("Read count after diamond subset: "+str(ds)+'\n')
        out.write("Eukaryota Read_count Protein_count Total_Protein_Seq_Len: "+", ".join(t)+'\n')


