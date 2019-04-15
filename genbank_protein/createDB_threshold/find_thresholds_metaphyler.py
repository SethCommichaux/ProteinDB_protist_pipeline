import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-nodes", help="NCBI nodes.dmp file")
parser.add_argument("-names", help="NCBI names.dmp file")
parser.add_argument("-f", help="Protist functions file")
args = parser.parse_args()

names = {i.strip().split('\t')[1].upper():i.strip().split('\t')[0] for i in open(args.names)}
nodes = {i.strip().split('\t')[0]:i.strip().split('\t')[1] for i in open(args.nodes)}
lineage = {i.strip().split('\t')[0]:i.strip().split('\t')[3] for i in open(args.f)}

print 'Built dictionaries...'

def get_taxa_rank(x,y,rank):
        x = map(lambda zz:zz.strip(),x.split(';'))
        y = map(lambda zz:zz.strip(),y.split(';'))
        for i in x:
                if nodes[names[i]] == rank:
                        x = i
                        break
        else: x = None

        for i in y:
                if nodes[names[i]] == rank:
                        y = i
                        break
        else: y = None

        return x,y


def find_cutoff(positives,negative):
        if negative == []: return min(positives)
        cutoff = 100000
        count = 1000000
        tmp = 0
        for i in set(positives):
                for j in positives:
                        if j < i: tmp += 1
                for k in negative:
                        if k > i: tmp += 1
                if tmp < count:
                        count = tmp
                        cutoff = i
                else: return cutoff
        return cutoff


def linear_regression(cuts):
        regressions = {}
        for i,j in cuts.items():
                x,y = [],[]
                for k in j:
                        x.append(float(k[1]))
                        y.append(float(k[0]))
                x,y = np.array(x),np.array(y)
                if len(x) <= 1: continue
                n = np.size(x)
                m_x, m_y = np.mean(x), np.mean(y)
                SS_xy = np.sum(y*x) - n*m_y*m_x
                SS_xx = np.sum(x*x) - n*m_x*m_x
                slope = SS_xy / SS_xx
                y_int = m_y - slope*m_x
                regressions[i] = [slope,y_int]
        return regressions



results = {}
subject_scores = {}

for f in ['30_5.pep.diamond','60_10.pep.diamond','90_15.pep.diamond','120_20.pep.diamond']:
        bd,nd,step = 'busco'+f,'negatives'+f,int(f.split('_')[0])
        for rank in ['phylum','class','order','family','genus','species']:
                scores = {}
                for i in open(bd):
                        tmp = i.strip().split('\t')
                        query,subject,aln,bitscore = '_'.join(tmp[0].split('_')[:2]),'_'.join(tmp[1].split('_')[:2]),float(tmp[3]),float(tmp[11])
                        if aln < step: continue
                        q_lin = lineage[query]
                        s_lin = lineage[subject]
                        q_lin,s_lin = get_taxa_rank(q_lin,s_lin,rank)
                        if None in [q_lin,s_lin]: continue
                        if q_lin == s_lin:
                                if subject not in scores:
                                        scores[subject] = [[bitscore],[],s_lin]
                                else:
                                        scores[subject][0].append(bitscore)


                if subject_scores == {}:
                        for i in open(nd):
                                tmp = i.strip().split('\t')
                                query,subject,aln,bitscore = '_'.join(tmp[0].split('_')[:2]),'_'.join(tmp[1].split('_')[:2]),float(tmp[3]),float(tmp[11])
                                if aln < step: continue
                                if subject not in subject_scores: subject_scores[subject] = [bitscore]
                                else: subject_scores[subject].append(bitscore)
                                if subject in scores: scores[subject][1].append(bitscore)
                else:
                        for k,v in subject_scores.items():
                                if k in scores: scores[k][1] = v


                for k,v in scores.items():
                        positives = sorted(v[0])
                        negative = sorted(v[1])
                        cutoff = find_cutoff(positives,negative)
                        if k not in results:
                                results[k] = {v[2]+'_'+rank:[(cutoff,step)]}
                        else:
                                if v[2]+'_'+rank not in results[k]:
                                        results[k][v[2]+'_'+rank] = [(cutoff,step)]
                                else:
                                        results[k][v[2]+'_'+rank].append((cutoff,step))


with open('tester2','w') as out:
        out.write('UniRef100\tSlope\tY_intercept\n')
        for k,v in results.items():
                regressions = linear_regression(v)
                if regressions == {}: continue
                out.write(k+'\t'+str(regressions)+'\n')
