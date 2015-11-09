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
	results = {}
	for w in sorted(distances, key=distances.get):
		if count <= 5:
			results[w] = str(distances[w])
			count += 1
		else:
			break
	return results
if __name__ == "__main__":

	print "Loading database..."
	with open(sys.argv[2]) as data_file:
		data = json.load(data_file)

	#Calculate tetranucleotide frequencies for query
	print "finding test virus-host pairs..."
	pairs = []
	phage_names = []
	with open(sys.argv[1]) as virus_host:
		for line in virus_host:
			try:
				name = line.split("\t")[4].split()[0]
				for key in data.keys():
					if name in key.split():
						phage_names.append([name,line.split("\t")[4]])
						pairs.append([name,line.split("\t")[0].strip()])
						break
			except:
				print "err"
	print pairs

	i = 0
	for pair in pairs:
		seqf = './viruses/' + pair[1]
		print seqf + "," + pair[1]
		handle = open(seqf, "rU")
		records = list(SeqIO.parse(handle, "fasta"))
		handle.close()
		seqs = records[0]
		print ""
		print "Calculating tetranucleotide frequencies for query..."
		tetramers = calc_tetra(seqs)

		#Compare
		
		print "Comparing query tetranucleotide frequencies to database..."
		print seqs.id + ":"
		hits = read_data(tetramers, data)
		found = False
		for hit in hits.keys():
			if pair[0] in hit:
				print "Success," + "," + str(phage_names[i]) + "," + str(pair) + "," + str(hits)
				found = True
				break
		if not found:
			print "Fail," + "," + str(phage_names[i]) + "," + str(pair) + "," + str(hits)

		i += 1