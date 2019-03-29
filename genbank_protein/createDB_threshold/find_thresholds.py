# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d", help="Path to protists.pep diamond index file")
parser.add_argument("-m", help="Path to protist_functions.txt file")
args = parser.parse_args()

def build_dicts():
        uniref100_2_Lineage = {}
        for i in open(args.m):
                tmp = i.strip().split('\t')
                uniref100 = tmp[0]
                lineage = tmp[3]
                uniref100_2_Lineage[uniref100] = lineage
        return uniref100_2_Lineage


def get_LCA(x,y):
        c=0
        for i in range(min(len(x),len(y))):
                if x[i] == y[i]:
                        c += 1
        new_lineage = x[:c]
        new_lineage = ';'.join(new_lineage)+';'
        return new_lineage


def lineage_thresholds(lineages,pidents,lineage,thresholds_V):
        line = ''
        threshold_by_lineage = []
        tmp = []
        for i in pidents:
                if i > 99.9: tmp.append(99.9)
                else: tmp.append(i)
        pidents = tmp
        for i in lineage.split(';')[:-1]:
                line += i+';'
                if line not in lineages:
                        threshold_by_lineage.append('NA')
                else:
                        threshold_by_lineage.append(min(pidents))
                        cutoff = max(thresholds_V[line])
                        tmp = [i for i in pidents if i > cutoff]
                        pidents = tmp
                        if pidents == []: pidents = [99]

        return threshold_by_lineage


def purge_results(results):
        thresholds = {}

        for k,v in results.items():
                lineage = uniref100_2_Lineage['_'.join(k.split('_')[:2])]
                lineage = ';'.join(map(lambda x:x.strip(),lineage.split(';')))
                if len(v) == 1:
                        thresholds[k] = {lineage+';':[99]}
                else:
                        lineage_chunks = map(lambda x: x.strip(), lineage.split(';'))
                        for z,w in v.items():
                                lineage2 = map(lambda x: x.strip(), z.split(';'))
                                LCA = get_LCA(lineage_chunks,lineage2)
                                if k not in thresholds:
                                        thresholds[k] = {LCA:sorted(w)}
                                else:
                                        if LCA not in thresholds[k]:
                                                thresholds[k][LCA] = sorted(w)
                                        else:
                                                thresholds[k][LCA] += w
                                                thresholds[k][LCA] = sorted(thresholds[k][LCA])

        with open('test2','a') as out:
                lineage = uniref100_2_Lineage['_'.join(k.split('_')[:2])]
                lineage = ';'.join(map(lambda x:x.strip(),lineage.split(';')))+';'
                for k,v in sorted(thresholds.items()):
                        pidents = []
                        lineages = []
                        for z,w in v.items():
                                pidents += w
                                lineages.append(z)
                        threshold_by_lineage = lineage_thresholds(lineages,pidents,lineage,v)
                        out.write('_'.join(k.split('_')[:2])+'\t'+k.split('_')[-1]+'\t'+lineage+'\t'+str(threshold_by_lineage)+'\n')



uniref100_2_Lineage = build_dicts()

results = {}
singletons = 0
switch = ''

for i in open(args.d):
        tmp = i.strip().split('\t')
        uniref100,hit,pident = tmp[0],tmp[1],float(tmp[2])
        id = '_'.join(uniref100.split('_')[:2])
        if id != switch:
                if switch == '':
                        switch = id
                else:
                        switch = id
                        purge_results(results)
                        results = {}
        lineage = ';'.join(map(lambda x:x.strip(),uniref100_2_Lineage[hit].split(';')))+';'
        if uniref100 not in results:
                results[uniref100] = {lineage:[pident]}
        else:
                if lineage not in results[uniref100]:
                        results[uniref100][lineage] = [pident]
                else:
                        results[uniref100][lineage].append(pident)


