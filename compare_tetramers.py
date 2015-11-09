import json
from scipy.spatial import distance
import sys
from Bio import SeqIO

def calc_tetra(seqs):

	tetramers = {}
	for a in ['A', 'C', 'G', 'T']:
		for b in ['A', 'C', 'G', 'T']:
			for c in ['A', 'C', 'G', 'T']:
				for d in ['A', 'C', 'G', 'T']:
					tetramers[a+b+c+d] = 0

	for seq in seqs:
		start = 0
		end = 4	
		for i in range(0,len(str(seqs.seq))):
			if len(str(seqs.seq[start:end])) == 4:
				try:
					tetramers[str(seqs.seq[start:end])] += 1
				except:
					pass	
			start += 1
			end += 1
	
	return tetramers


def read_data(tetramers1, tetramers2):
	#Normalize

	total = sum(tetramers1.values())
	for k in tetramers1.keys():
		tetramers1[k] = float(tetramers1[k]) / float(total)

	total = sum(tetramers2.values())
	for k in tetramers2.keys():
		tetramers2[k] = float(tetramers2[k]) / float(total)

	query_dat = []
	for d in sorted(tetramers1.keys()):
		query_dat.append(tetramers1[d])

	subject_dat = []
	for d in sorted(tetramers2.keys()):
		subject_dat.append(tetramers2[d])

	
	print "Calculated euclidean distance: " + str(round(distance.euclidean(query_dat, subject_data),5))

if __name__ == "__main__":

	#Calculate tetranucleotide frequencies for query
	
	handle = open(sys.argv[1], "rU")
	records = list(SeqIO.parse(handle, "fasta"))
	handle.close()

	handle = open(sys.argv[2], "rU")
	records2 = list(SeqIO.parse(handle, "fasta"))
	handle.close()

	print ""
	print "Calculating tetranucleotide frequencies for query..."
	tetramers = calc_tetra(records)
	tetramers2 = calc_tetra(records2)
	#Compare
	
	print "Comparing query and subject tetranucleotide frequencies..."
	read_data(tetramers, tetramers2)