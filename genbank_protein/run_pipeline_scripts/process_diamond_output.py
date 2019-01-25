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
read_LCA,species_protein_count,bitscores = {},{},{}


# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen
for i in open(args.d):
	tmp = i.strip().split('\t')
	read,gb_id,pident,qlen,start,stop,bitscore = tmp[0],tmp[1],float(tmp[2]),int(tmp[-1]),int(tmp[6]),int(tmp[7]),float(tmp[11])
	if gb_id not in fullnamelineage: continue
	lineage = fullnamelineage[gb_id]
	if lineage in [None,'']: continue
	if 'UNCLASSIFIED SEQUENCES' in lineage: continue
	if 'CELLULAR ORGANISMS' not in lineage: continue
	if read not in bitscores: 
		bitscores[read] = bitscore
	else:
		if bitscore < bitscores[read]: continue
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


taxaGroup = Counter(read_LCA.values())

for k,v in species_protein_count.items():
	if k in taxaGroup:
		if len(v) < 5:
			tmp = taxaGroup[k]
			tmp2 = ';'.join(k.split(';')[:-2])+';'
			del taxaGroup[k]
			taxaGroup[tmp2] = tmp

with open(args.d+'.results','w') as out:
	for k,v in sorted(taxaGroup.items(),key=lambda x: x[1]):
		if k == None: continue
		if 'EUKARYOTA' in k:
			if 'METAZOA' not in k:
				if 'FUNGI;' not in k:
					if 'EMBRYOPHYTA' not in k:
						out.write(str(k)+'\t'+str(v)+'\t'+str(float(v)/read_count)+'\n')

