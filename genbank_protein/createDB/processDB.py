import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-g", help="Path to genbank protein fasta file")
parser.add_argument("-t", help="path to NCBI fullnamelineage.dmp taxonomy file")
args = parser.parse_args()

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


taxaName2_ID_Lineage = {}

for i in open(args.t):
        tmp = i.strip().split('\t|\t')
        taxID = tmp[0]
        taxaName = tmp[1].upper().strip('\t|\t')
        lineage = tmp[2].upper().strip('\t|\t')
        taxaName2_ID_Lineage[taxaName] = [taxID,lineage+taxaName]

print "Parsed NCBI taxa lineage file!!!"


with open('binningDB.pep','w') as out, open('genbank_map_ncbi_taxonomy.txt','w') as out2, open('queryDB.pep','w') as out3, open('environmental_eukaryote.fasta','w') as out4:
        out2.write('GENBANKID\tTAXANAME\tTAXAID\tLINEAGE\tPROTEIN_NAME\tPROTEIN_LENGTH\tFUNCTION\n')
        for i in SeqIO.parse(args.g,'fasta'):
                d = str(i.description)
                id = str(i.id)
                s = str(i.seq)
                l = len(s)
                try:
                        taxa = d.split('[')[1].split(']')[0].upper()
                        lineage = taxaName2_ID_Lineage[taxa][1]
                        if "CELLULAR ORGANISMS" not in lineage:
                                if 'VIRUS' not in lineage: continue
                        if 'EUKARYOTA' in lineage:
                                if 'ENVIRONMENTAL SAMPLES' in lineage:
                                        out4.write(">"+d+"\n"+s+"\n")
                                        continue
                                if 'UNCLASSIFIED EUKARYOTES' in lineage:
                                        out4.write(">"+d+"\n"+s+"\n")
                                        continue
                        taxaID = taxaName2_ID_Lineage[taxa][0]
                        proteinName = ' '.join(d.split(' ')[1:]).split('[')[0]
                        out2.write(id+'\t'+taxa+'\t'+taxaID+'\t'+lineage+'\t'+proteinName+'\t'+l+'\t'+'NA'+'\n')
                        r = isProtist(taxa)
                        if r == True:
                                out.write(">"+d+"\n"+s+"\n")
                        out3.write(">"+d+"\n"+s+"\n")
                except: continue
