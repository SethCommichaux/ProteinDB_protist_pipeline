import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-g", help="path to input genbank file of proteins: fasta format")
parser.add_argument("-t", help="path to NCBI fullnamelineage.dmp taxonomy file")
parser.add_argument("-o", help="desired output file name")
args = parser.parse_args()

taxaName2_ID_Lineage = {}

for i in open(args.t):
	tmp = i.strip().split('\t|\t')
	taxID = tmp[0]
	taxaName = tmp[1].upper().strip('\t|\t')
	lineage = tmp[2].upper().strip('\t|\t')
	taxaName2_ID_Lineage[taxaName] = [taxID,lineage+taxaName]

with open(args.o,'w') as out:
	out.write('GENBANKID\tTAXANAME\tTAXAID\tLINEAGE\tPROTEIN_NAME\n')
	for i in SeqIO.parse(args.g,'fasta'):
		d = str(i.description)
		id = str(i.id)
		s = str(i.seq)
		if '[' not in d:
			out.write(id+'\t'+taxa+'\t'+'\t'+'\n')
			continue
		proteinName = ' '.join(d.split(' ')[1:]).split('[')[0]
		taxa = d.split('[')[1].split(']')[0].upper()
		if taxa not in taxaName2_ID_Lineage:
			out.write(id+'\t'+taxa+'\t'+'\t'+'\t'+proteinName+'\n')
		else:
			out.write(id+'\t'+taxa+'\t'+taxaName2_ID_Lineage[taxa][0]+'\t'+taxaName2_ID_Lineage[taxa][1]+'\t'+proteinName+'\n')
