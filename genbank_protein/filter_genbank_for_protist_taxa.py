import argparse
from Bio import SeqIO

def isProtist(taxa):
	result = False
	if taxa in lineage:
		l = lineage[taxa]
		if 'EUKARYOTA' in l:
			if 'EMBRYOPHYTA' not in l:
				if 'METAZOA' not in l:
					if 'FUNGI' not in l:
						result = True 
	return result	


parser = argparse.ArgumentParser()
parser.add_argument("-g", help="Path to genbank protein fasta file")
parser.add_argument("-l", help="Path to NCBI fullnamelineage.dmp file")
parser.add_argument("-o", help="Path name for output fasta file of protist proteins")
args = parser.parse_args()


lineage = {i.strip().split('\t|\t')[1].upper():i.strip().split('\t|\t')[2].upper() for i in open(args.l)}
print "Parsed NCBI taxa lineage file!!!"


with open(args.o,'w') as out, open(args.g+'without_protists','w') as out2:
	for i in SeqIO.parse(args.g,'fasta'):
		d = str(i.description)
		try:
			taxa = d.split('[')[1].split(']')[0].upper()
			s = str(i.seq)
			r = isProtist(taxa)
			if r == False:
				out2.write(">"+d+"\n"+s+"\n") 
			else:
				out.write(">"+d+"\n"+s+"\n")
		except:
			continue

