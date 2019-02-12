import argparse
import sys
from collections import Counter

parser = argparse.ArgumentParser()
parser.add_argument("-d", type=str,
						help="Diamond output file (--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen)")
parser.add_argument("-f", type=str,
						help="Fasta/Fastq reads file. Suffix must be either *.fasta or *.fastq")
parser.add_argument("-t", type=str,
						help="The genbank_map_ncbi_taxonomy.txt file from creating database.")
parser.add_argument("-p", type=str,
                                                help="Text file that maps protein IDs to protein lengths in amino acids")
args = parser.parse_args()

def count_reads(x):
	c = 0
	if x.endswith('.fasta'):
		for i in open(x):
			if i[0] == ">": c += 1
		return c
	elif x.endswith('.fastq'):
		for i in open(x):
			if i[0] == "+": c += 1
		return c
	else:
		print "Process Terminated: Need *.fasta or *.fastq file!!!"
		sys.exit()


def build_taxa_dicts():
	fullnamelineage,gbID_species = {},{}
	for i in open(args.t):
		tmp = i.strip().split('\t')
		try:
			gb_id = tmp[0]
			taxaName = tmp[1].upper()
			taxaID = tmp[2]
			lineage = tmp[3]
			proteinName = tmp[4]
			fullnamelineage[gb_id] = lineage
			gbID_species[gb_id] = taxaName
		except: continue
	return fullnamelineage,gbID_species


def build_protein_len_dict():
	proteinID_len = {i.strip().split('\t')[0]:int(i.strip().split('\t')[1]) for i in open(args.p)}
	return proteinID_len


def LCA(x,y):
	new_lineage = ''
	x=map(lambda z: z.strip(), x.split(';'))
	y=map(lambda z: z.strip(), y.split(';'))
	for i in range(min(len(x),len(y))):
		if x[i] == y[i]:
			new_lineage += x[i]+';'
		else:
			return new_lineage


read_count = count_reads(args.f)
fullnamelineage,gbID_species = build_taxa_dicts()
proteinID_len = build_protein_len_dict()
read_LCA,species_protein_count = {},{}


# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen
for i in open(args.d):
	tmp = i.strip().split('\t')
	read,gb_id,pident,qlen,start,stop,bitscore = tmp[0],tmp[1],float(tmp[2]),int(tmp[-1]),int(tmp[6]),int(tmp[7]),float(tmp[11])
	if gb_id not in fullnamelineage: continue
	lineage = fullnamelineage[gb_id]
	if 'CELLULAR ORGANISMS' not in lineage: continue
	if 'EUKARYOTA' in lineage:
		if 'ENVIRONMENTAL SAMPLES' in lineage: continue
	if read not in read_LCA:
		read_LCA[read] = lineage
	else:
		if read_LCA[read] == None: continue
		read_LCA[read] = LCA(read_LCA[read],lineage)
	if lineage not in species_protein_count:
		species_protein_count[lineage] = [gb_id]
	else:
		if gb_id not in species_protein_count[lineage]:
			species_protein_count[lineage].append(gb_id)

fullnamelineage = {}

species_protein_count_cumulative = {}

for species,protein_count in species_protein_count.items():
	tmp = map(lambda z: z.strip(),species.split(';'))
	for rank in range(1,len(tmp)):
		line = ';'.join(tmp[0:rank])+';'
		if line in species_protein_count_cumulative:
			species_protein_count_cumulative[line] += protein_count
		else:
			species_protein_count_cumulative[line] = protein_count

for k,v in species_protein_count_cumulative.items():
	species_protein_count_cumulative[k] = set(v)


species_protein_count = {}

LCA_readcount = Counter(read_LCA.values())
read_LCA = {}
read_LCA_cumulative = {}

for LCA,read_count in LCA_readcount.items():
	if LCA == None: continue
	tmp = map(lambda z: z.strip(),LCA.split(';'))
	for rank in range(1,len(tmp)):
		line = ';'.join(tmp[0:rank])+';'
		if line in read_LCA_cumulative:
			read_LCA_cumulative[line] += read_count
		else:
			read_LCA_cumulative[line] = read_count


LCA_readcount = {}

with open(args.d+'.results','w') as out:
	for k,v in sorted(read_LCA_cumulative.items(),key=lambda x: x[1]):
		if k == None: continue
		if v < 5: continue
		if len(species_protein_count_cumulative[k]) < 5: continue
		if 'EUKARYOTA' in k:
			if 'METAZOA' not in k:
				if 'FUNGI;' not in k:
					if 'EMBRYOPHYTA' not in k:
						out.write(str(k)+'\t'+str(v)+'\t'+str(float(v)/read_count)+'\t'+str(len(species_protein_count_cumulative[k]))+'\t'+str(sum(map(lambda z: proteinID_len[z],species_protein_count_cumulative[k])))+'\n')

