import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-g", help="Path to uniprot protein fasta file")
parser.add_argument("-t", help="path to NCBI fullnamelineage.dmp taxonomy file")
args = parser.parse_args()


# parse NCBI taxonomy fullnamelineage.dmp file for taxaIDs, taxaNames and taxonomic lineages
def build_taxa_dict(): 
	taxaID2_Name_Lineage = {}
	for i in open(args.t):
		tmp = i.strip().split('\t|\t')
		taxaID = tmp[0]
		taxaName = tmp[1].upper().strip('\t|\t')
		lineage = tmp[2].upper().strip('\t|\t')
		taxaID2_Name_Lineage[taxaID] = [taxaName,lineage+taxaName]
	print "Taxa dictionary built!!!"
	return taxaID2_Name_Lineage

taxaID2_Name_Lineage = build_taxa_dict()


with open('binningDB.pep','w') as out, open('nonprotists.pep','w') as out2:
	for i in open(args.g):
		if i.startswith(">"):
			query = False
			UniRef100 = i.strip().split(' ')[0][1:]
			taxaID = i.split('TaxID=')[1].split(' ')[0]
			if taxaID not in taxaID2_Name_Lineage: continue
			lineage = taxaID2_Name_Lineage[taxaID][1]
			if 'BACTERIA' in lineage: out2.write(i)
			elif 'VIRUS' in lineage: out2.write(i)
			elif 'ARCHAEA' in lineage: out2.write(i)
			elif 'EUKARYOTA' in lineage: 
				if 'EMBRYOPHYTA' in lineage: out2.write(i)
				elif 'METAZOA' in lineage: out2.write(i)
				elif 'FUNGI' in lineage: out2.write(i)
				elif 'ENVIRONMENTAL SAMPLES' in lineage: continue
				elif 'UNCLASSIFIED EUKARYOTES' in lineage: continue
				else:
					out.write(i)
					query = True
		elif query == False: out2.write(i)
		else: out.write(i)
