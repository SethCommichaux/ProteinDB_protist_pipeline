import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-k", help="Kaiju output file")
parser.add_argument("-np", help="Non-protists fasta file")
args = parser.parse_args()


kaiju = {i.strip().split('\t')[1].split(' ')[0]:0 for i in open(args.k)}

with open('queryDB.pep','a') as out:
	for i in open(args.np):
		if i.startswith(">"):
			switch = False
			id = i.split(' ')[0][1:]
			if id in kaiju:
				out.write(i)
				switch = True
		elif switch == True: out.write(i)
