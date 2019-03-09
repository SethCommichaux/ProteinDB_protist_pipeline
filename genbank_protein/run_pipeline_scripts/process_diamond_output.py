import argparse
import sys
from collections import Counter

parser = argparse.ArgumentParser()
parser.add_argument("-d", type=str,
						help="Diamond output file (--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen)")
parser.add_argument("-t", type=str,
						help="The genbank_map_ncbi_taxonomy.txt file from creating database.")
parser.add_argument("-f", type=int,
			help="Read count of input file")
parser.add_argument("-num_reads", type=str,
						help="Minimum number of mapped reads to include phylogroup")
parser.add_argument("-num_prots", type=str,
						help="Minimum number of mapped proteins to include phylogroup")
args = parser.parse_args()

print "Parsed args!!!"

def build_taxa_dicts():
	fullnamelineage = {i.strip().split('\t')[1]:0 for i in open(args.d)}
	for i in open(args.t):
		tmp = i.strip().split('\t')
		gb_id = tmp[0]
		prot_len = int(tmp[5])
		lineage = tmp[3]
		fullnamelineage[gb_id] = [prot_len,lineage]
	print "Built taxa dictionary!!!"
	return fullnamelineage


def LCA(x,y):
	new_lineage = ''
	x=map(lambda z: z.strip(), x.split(';'))
	y=map(lambda z: z.strip(), y.split(';'))
	for i in range(min(len(x),len(y))):
		if x[i] == y[i]:
			new_lineage += x[i]+';'
		else:
			return new_lineage


fullnamelineage = build_taxa_dicts()
read_LCA = {}
species_protein_count = {}
ANI = {}
ANI_lineage = {}
ANI_cumulative = {}
species_protein_count_cumulative = {}
read_LCA_cumulative = {}


# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen
for i in open(args.d):
	tmp = i.strip().split('\t')
	read,gb_id,pident,qlen,start,stop,bitscore = tmp[0],tmp[1],float(tmp[2]),int(tmp[-1]),int(tmp[6]),int(tmp[7]),float(tmp[11])
	lineage = fullnamelineage[gb_id][1]
 	if 'VIRUS' in lineage: 
		read_LCA[read] = None
		continue
	if 'BACTERIA' in lineage: 
		read_LCA[read] = None
		continue
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
	if read not in ANI: ANI[read] = [pident]
	else: ANI[read].append(pident)

print "Processed diamond file!!!"


for read,pident in ANI.items():
	if read in read_LCA:
		if read_LCA[read] in ANI_lineage:
			ANI_lineage[read_LCA[read]] += pident
		else:
			ANI_lineage[read_LCA[read]] = pident

ANI = {}

for lineage,pident in ANI_lineage.items():
	if lineage == None: continue
	tmp = map(lambda z: z.strip(),lineage.split(';'))
	for rank in range(1,len(tmp)):
		line = ';'.join(tmp[0:rank])+';'
		if line in ANI_cumulative:
			ANI_cumulative[line] += pident
		else:
			ANI_cumulative[line] = pident

ANI_lineage = {}

for species,protein_count in species_protein_count.items():
	tmp = map(lambda z: z.strip(),species.split(';'))
	for rank in range(1,len(tmp)):
		line = ';'.join(tmp[0:rank])+';'
		if line in species_protein_count_cumulative:
			species_protein_count_cumulative[line] = species_protein_count_cumulative[line] | set(protein_count)
		else:
			species_protein_count_cumulative[line] = set(protein_count)


LCA_readcount = Counter(read_LCA.values())
species_protein_count,read_LCA = {},{}


for LCA,read_count in LCA_readcount.items():
	if LCA == None: continue
	tmp = map(lambda z: z.strip(),LCA.split(';'))
	for rank in range(1,len(tmp)):
		line = ';'.join(tmp[0:rank])+';'
		if line in read_LCA_cumulative:
			read_LCA_cumulative[line] += read_count
		else:
			read_LCA_cumulative[line] = read_count

print "Proccessed dictionaries!!!"
LCA_readcount = {}

with open(args.d+'.results','w') as out:
	out.write('Phylogroup\tRead_count\tRel_abund\tProtein_count\tProt_Seq_Len\tANI\n')
	for k,v in sorted(read_LCA_cumulative.items(),key=lambda x: x[1]):
		if k == None: continue
		if v < int(args.num_reads): continue
		if len(species_protein_count_cumulative[k]) < int(args.num_prots): continue
		if 'EUKARYOTA' in k:
			if 'METAZOA' not in k:
				if 'FUNGI;' not in k:
					if 'EMBRYOPHYTA' not in k:
						out.write(str(k)+'\t'+str(v)+'\t'+str(v/float(args.f))+'\t'+str(len(species_protein_count_cumulative[k]))+'\t'+str(sum(map(lambda z: fullnamelineage[z][0],species_protein_count_cumulative[k])))+'\t'+str(float(sum(ANI_cumulative.get(k,1)))/len(ANI_cumulative.get(k,1)))+'\n')
