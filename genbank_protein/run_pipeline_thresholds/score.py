import sys
import ast
import argparse
import collections

parser = argparse.ArgumentParser()
parser.add_argument("-t", help="Protist protein-specific thresholds file")
parser.add_argument("-f", help="NCBI fullnamelineage.dmp taxonomy file")
parser.add_argument("-d", help="Diamond output file to be analyzed")
parser.add_argument("-o", help="Desired output file")
args = parser.parse_args()

lineages = {i.strip().split('\t|\t')[1].upper():i.strip().split('\t|\t')[2].strip('\t|').upper()+i.strip().split('\t|\t')[1].upper() for i in open(args.f)}

protein_thresholds = args.t
diamond_file = args.d

thresholds = {}

for h,i in enumerate(open(protein_thresholds)):
        if h == 0: continue
        tmp = i.strip().split('\t')
        uniref100 = tmp[0]
        cutoffs = ast.literal_eval(tmp[1])
        thresholds[uniref100] = cutoffs


results = {}
multiprotist_reads = collections.Counter([i.strip().split('\t')[0] for i in open(diamond_file)])

results_cumulative = {}

with open(args.o,'w') as out:
        for i in open(diamond_file):
                result = {}
                best = 0
                tmp = i.strip().split('\t')
                read = tmp[0]
                uniref100 = tmp[1]
                aln = float(tmp[3])
                bitscore = float(tmp[11])
                multi = multiprotist_reads[read]
                if multi > 1: uniq = 0
                else: uniq = 1
                if uniref100 in thresholds:
                        for k,v in thresholds[uniref100].items():
                                if bitscore >= (v[1] + v[0]*aln):
                                        result[k.split('_')[1]] = k.split('_')[0]
                        if result == {}: continue
                        else:
                                for z in ['species','genus','family','order','class','phylum']:
                                        if z in result:
                                                taxa = result[z]+'_'+z
                                                if taxa in results:
                                                        results[taxa][0] += 1
                                                        results[taxa][1] += uniq
                                                else: results[taxa] = [1,uniq]
                                                break

        for k,v in sorted(results.items(),key=lambda x:x[1]):
                if v[1] > 0:
                        taxa = k.split('_')[0]
                        lineage = lineages.get(taxa,0).strip()
                        if lineage == 0:
                                continue
                        line = ''
                        for i in lineage.split(';'):
                                line += i+';'
                                if line not in results_cumulative:
                                        results_cumulative[line] = v[0]
                                else:
                                        results_cumulative[line] += v[0]

        for k,v in sorted(results_cumulative.items()):
                out.write(k+'\t'+str(v)+'\n')

