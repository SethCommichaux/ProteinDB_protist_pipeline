import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-g", help="Path to genbank protein fasta file")
parser.add_argument("-t", help="path to NCBI fullnamelineage.dmp taxonomy file")
args = parser.parse_args()

euk,protists,enviro_euk,unclassified_euk,genbank,bacteria,virus,archaea,plant,fungi,animal = 0,0,0,0,0,0,0,0,0,0,0
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
		genbank += 1
                d = str(i.description)
                id = str(i.id)
                s = str(i.seq)
                l = str(len(s))
                taxa = d.split('[')[-1].split(']')[0].upper()
		if taxa not in taxaName2_ID_Lineage: continue
                lineage = taxaName2_ID_Lineage[taxa][1]
                taxaID = taxaName2_ID_Lineage[taxa][0]
                proteinName = ' '.join(d.split(' ')[1:]).split('[')[0]
                out2.write(id+'\t'+taxa+'\t'+taxaID+'\t'+lineage+'\t'+proteinName+'\t'+l+'\t'+'NA'+'\n')
                out3.write(">"+d+"\n"+s+"\n")
                if 'BACTERIA' in lineage: bacteria +=1
                elif 'VIRUS' in lineage: virus +=1
                elif 'ARCHAEA' in lineage: archaea +=1
                elif 'EUKARYOTA' in lineage: 
                    euk +=1
                    if 'EMBRYOPHYTA' in lineage: plant +=1
                    elif 'METAZOA' in lineage: animal +=1
                    elif 'FUNGI' in lineage: fungi +=1
                    elif 'ENVIRONMENTAL SAMPLES' in lineage: enviro_euk +=1; out4.write(">"+d+"\n"+s+"\n")
                    elif 'UNCLASSIFIED EUKARYOTES' in lineage: unclassified_euk +=1; out4.write(">"+d+"\n"+s+"\n")
                    else: protists +=1; out.write(">"+d+"\n"+s+"\n")

with open('Genbank_summary.txt','w') as out5:
	out5.write('Total: '+str(genbank)+'\n'+'Virus: '+str(virus)+'\n'+'Bacteria: '+str(bacteria)+'\n'+'Archaea: '+str(archaea)+'\n'+'Eukaryote: '+str(euk)+'\n'+'Plant: '+str(plant)+'\n'+'Fungi: '+str(fungi)+'\n'+'Animal: '+str(animal)+'\n'+'Protists: '+str(protists)+'\n'+'Environmental_eukaryotes: '+str(enviro_euk)+'\n'+'Unclassified_eukaryotes: '+str(unclassified_euk)+'\n')
