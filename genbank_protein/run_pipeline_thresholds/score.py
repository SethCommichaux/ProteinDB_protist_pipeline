import sys
import ast
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t", help="Protist protein-specific thresholds file")
parser.add_argument("-d", help="Diamond output file to be analyzed")
parser.add_argument("-o", help="Desired output file")
args = parser.parse_args()

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
#multiprotist_reads = {}

with open(args.o,'w') as out:
        for i in open(diamond_file):
                result = {}
                best = 0
                tmp = i.strip().split('\t')
                read = tmp[0]
                uniref100 = tmp[1]
                aln = float(tmp[3])
                bitscore = float(tmp[11])
                if uniref100 in thresholds:
                        for k,v in thresholds[uniref100].items():
                                if bitscore >= (v[1] + v[0]*aln):
                                        result[k.split('_')[1]] = k.split('_')[0]
                        if result == {}: continue
                        else:
                                for z in ['species','genus','family','order','class','phylum']:
                                        if z in result:
                                                taxa = result[z]+'_'+z
#                                               if read not in multiprotist_reads:
#                                                       multiprotist_reads[read] = [taxa]
#                                               else: multiprotist_reads[read].append(taxa)
                                                if taxa in results:
                                                        results[taxa] += 1
                                                else: results[taxa] = 1
                                                break
        for k,v in sorted(results.items(),key=lambda x:x[1]):
                out.write(k+'\t'+str(v)+'\n')


#for k,v in multiprotist_reads.items():
#       v=list(set(v))
#       if len(v) > 1: print k,'\t',v
