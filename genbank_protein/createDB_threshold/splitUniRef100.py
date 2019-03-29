import argparse
from Bio import SeqIO


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
                taxaID2_Name_Lineage[taxaID] = lineage+taxaName
        print "Taxa dictionary built!!!"
        return taxaID2_Name_Lineage

taxaID2_Name_Lineage = build_taxa_dict()


with open('protists.pep','w') as out:
        for i in SeqIO.parse(args.g,'fasta'):
                d = str(i.description)
                s = str(i.seq)
                taxaID = d.split('TaxID=')[1].split(' ')[0]
                if taxaID not in taxaID2_Name_Lineage: continue
                lineage = taxaID2_Name_Lineage[taxaID]
                if 'EUKARYOTA' in lineage:
                        if 'EMBRYOPHYTA' not in lineage:
                                if 'METAZOA' not in lineage:
                                        if 'FUNGI' not in lineage:
                                                if 'ENVIRONMENTAL SAMPLES' not in lineage:
                                                        if 'UNCLASSIFIED EUKARYOTES' not in lineage:
                                                                out.write(">"+d+"\n"+s+"\n")
