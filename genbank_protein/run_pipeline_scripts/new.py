import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d", type=str,
						help="Diamond output file (--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen)")
parser.add_argument("-t", type=str,
						help="The genbank_map_ncbi_taxonomy.txt file from creating database.")
parser.add_argument("-f", type=int,
						help="Read count of input file")
parser.add_argument("-ANI", type=float,
						help="Minimum average nucleotide identity (e.g. 80 for 80%) for reads aligned to phylogroup")
parser.add_argument("-num_reads", type=str,
						help="Minimum number of mapped reads to include phylogroup")
parser.add_argument("-num_prots", type=str,
						help="Minimum number of mapped proteins to include phylogroup")
args = parser.parse_args()


def build_taxa_dicts():
	fullnamelineage = {i.strip().split('\t')[1]:0 for i in open(args.d)}
	for h,i in enumerate(open(args.t)):
		if h == 0: continue
		tmp = i.strip().split('\t')
		gb_id = tmp[0]
		prot_len = int(tmp[5])
		lineage = tmp[3]
		fullnamelineage[gb_id] = [prot_len,lineage]
	print "Built taxa dictionary!!!"
	return fullnamelineage


def LCA(taxa):
	if len(taxa) == 1:
		return taxa[0]
	for i,j in enumerate(taxa):
		if 'EUKARYOTA' not in j:
			del first_pass[read]
			return False
		else:
			c = 0
			j = map(lambda z: z.strip(), x.split(';'))
			if i == 0:
				new_lineage = j
			else:
				for i in range(min(len(j),len(new_lineage))):
					if j[i] == new_lineage[i]:
						c += 1
				new_lineage = new_lineage[:c]
	new_lineage = ';'.join(new_lineage)+';'
	return new_lineage


def calc_ANI(pident,aln):
	total_len = sum(aln)
	weight = 0
	for i in range(len(pident)):
		weight += pident[i]*aln[i]/100.0
	if weight/total_len < args.ANI: return False
	return weight/total_len


fullnamelineage = build_taxa_dicts()
first_pass = {}

# collect values in file
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen
for i in open(args.d):
	tmp = i.strip().split('\t')
	read,gb_id,pident,aln = tmp[0],tmp[1],float(tmp[2]),float(tmp[3])
	lineage = fullnamelineage[gb_id][1]
	if read not in first_pass:
		first_pass[read] = [[lineage],[gb_id],[pident],[aln]]
	else:
		first_pass[read][0].append(lineage)
		first_pass[read][1].append(gb_id)
		first_pass[read][2].append(pident)
		first_pass[read][3].append(aln)


# calculate LCA for each read
for read,v in first_pass.items():
	taxa = list(set(v[0]))
	read_LCA = LCA(taxa)
	if read_LCA == False: 
		del first_pass[read]
		continue
	else: first_pass[read][0] = read_LCA


# group results by phylogroup
results = {}

for read,v in first_pass.items():
	LCA = v[0]
	proteins = v[1]
	pident = v[2]
	aln = v[3]
	if LCA not in results:
		results[LCA] = [1,proteins,pident,aln]
	else:
		results[LCA][0] += 1 # read_count
		results[LCA][1] += proteins
		results[LCA][2] += pident
		results[LCA][3] += aln


# calculate ANI for phylogroup
first_pass = {}

for LCA,v in results.items():
	read_count = v[0]
	protein_count = set(v[1])
	pident = v[2]
	aln = v[3]
	ANI = calc_ANI(pident,aln)
	if ANI == False:
		del results[LCA]
		continue
	results[LCA] = [read_count,protein_count,ANI,pident,aln]


# retrieve cumulative results for each phylogroup
cumulative_results = {}

for lineage,v in results.items():
	tmp = map(lambda z: z.strip(),lineage.split(';'))
	for rank in range(1,len(tmp)):
		line = ';'.join(tmp[0:rank])+';'
		if line not in cumulative_results:
			cumulative_results[line] = v
		else:
			cumulative_results[line][0] += v[0]
			cumulative_results[line][1] = cumulative_results[line][1] | v[1]
			pident += cumulative_results[line][3]
			aln += cumulative_results[line][4]
			cumulative_results[line][2] = calc_ANI(pident,aln)


with open(args.d+'.results','w') as out, open(args.d+'.functions','w') as out2:
	out.write('Phylogroup\tRead_count\tRel_abund\tProtein_count\tProt_Seq_Len\tANI\n')
	for k,v in sorted(cumulative_results.items()):
		if v[1] < int(args.num_reads): continue
		if len(v[2]) < int(args.num_prots): continue
		if 'EUKARYOTA' in k:
			if 'METAZOA' not in k:
				if 'FUNGI;' not in k:
					if 'EMBRYOPHYTA' not in k:
						out.write(str(k)+'\t'+str(v[0])+'\t'+str(v[0]/float(args.f))+'\t'+str(len(v[1]))+'\t'+str(sum(map(lambda z: fullnamelineage[z][0],v[1])))+'\t'+str(v[2])+'\n')
						out2.write(str(k)+'\t'+str(v[1])+'\n')
