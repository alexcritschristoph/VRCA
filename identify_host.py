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


def read_data(tetramers, data):
	#Normalize

	total = sum(tetramers.values())
	for k in tetramers.keys():
		tetramers[k] = float(tetramers[k]) / float(total)

	for species in data.keys():
		total = sum(data[species].values())
		for k in data[species].keys():
			data[species][k] = float(data[species][k]) / float(total)

	#compare	
	query_dat = []
	for d in sorted(tetramers.keys()):
		query_dat.append(tetramers[d])

	distances = {}
	for species in data.keys():
		subject_data = []
		for d in sorted(data[species].keys()):
			subject_data.append(data[species][d])
		distances[species] = round(distance.euclidean(query_dat, subject_data),5)
	
	count = 0
	for w in sorted(distances, key=distances.get):
		if count <= 5:
			print w + ":" + str(distances[w])
			count += 1
		else:
			break

if __name__ == "__main__":

	#Calculate tetranucleotide frequencies for query
	
	handle = open(sys.argv[1], "rU")
	records = list(SeqIO.parse(handle, "fasta"))
	handle.close()
	print "Loading database..."
	with open(sys.argv[2]) as data_file:
		data = json.load(data_file)

	for seqs in records:
		print ""
		print "Calculating tetranucleotide frequencies for query..."
		tetramers = calc_tetra(seqs)

		#Compare
		
		print "Comparing query tetranucleotide frequencies to database..."
		print seqs.id + ":"
		read_data(tetramers, data)