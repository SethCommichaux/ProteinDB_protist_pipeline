results = {}

for i in open('kaiju.fasta.diamond.subset'):
	tmp = i.strip().split('\t')
	protein = tmp[1]
	pident = float(tmp[2])
	sstart = float(tmp[8])
	if protein not in results:
		results[protein] = [[20*round(sstart/20.)],[pident]]
	else:
		results[protein][0].append(20*round(sstart/20.))
		results[protein][1].append(pident)


output = {}

with open('practice','w') as out:
	for j in open('/lustre/projects/SethCommichaux/ProteinDB_protist_pipeline/protist_thresholds.txt'):
		tmp = j.strip().split('\t')
		protein = tmp[0]
		start = int(tmp[1])
		lineage = tmp[2]
		thresholds = tmp[3].replace('[','').replace(']','').split(',')
		pick = ''
		if protein in results:
			if start in results[protein][0]:
				pident = results[protein][1][results[protein][0].index(start)]
				for h,i in enumerate(thresholds):
					if i.strip() != "'NA'":
						if pident > float(i):
							pick = h+1
				if pick == '': continue
				else:
					answer = ';'.join(lineage.split(';')[:pick])
					if answer in output: output[answer] += 1
					else: output[answer] = 1
	
	for k,v in sorted(output.items(),key=lambda x: x[1]):
		out.write(k+'\t'+str(v)+'\n')


