import sys
import ast
import argparse
import collections

def LCA(taxa):
        for i,j in enumerate(taxa):
                c = 0
                j = j.split(';')
                if i == 0:
                        new_lineage = j
                else:
                        for i in range(min(len(j),len(new_lineage))):
                                if j[i] == new_lineage[i]:
                                        c += 1
                        new_lineage = new_lineage[:c]
        new_lineage = ';'.join(new_lineage)
        return new_lineage


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


results, LCA_results, results_cumulative = {},[],{}

with open(args.o,'w') as out:
        for i in open(diamond_file):
                result = {}
                best = 0
                tmp = i.strip().split('\t')
                read = tmp[0]
                uniref100 = tmp[1]
                aln = float(tmp[3])
                if aln < 30: continue
                bitscore = float(tmp[11])
                uniq = 1 # need to remove line
                if uniref100 in thresholds:
                        for k,v in thresholds[uniref100].items():
                                if bitscore >= (v[1] + v[0]*aln):
                                        result[k.split('_')[1]] = k.split('_')[0]
                        if result == {}: continue
                        else:
                                for z in ['species','genus','family','order','class','phylum']:
                                        if z in result:
                                                taxa = result[z]
                                                lineage = map(lambda z:z.strip(),lineages[taxa].split(';'))
                                                lineage = ';'.join(lineage)
                                                if read not in results:
                                                        results[read] = [lineage]
                                                else: results[read].append(lineage)
                                                break

        for k,v in results.items():
                if len(v) == 1:
                        LCA_results.append(v[0])
                else: LCA_results.append(LCA(v))

        for k,v in collections.Counter(LCA_results).items():
                line = ''
                for i in k.split(';'):
                        line += i+';'
                        if line not in results_cumulative:
                                results_cumulative[line] = v
                        else:
                                results_cumulative[line] += v

        for k,v in sorted(results_cumulative.items()):
                out.write(k+'\t'+str(v)+'\n')


