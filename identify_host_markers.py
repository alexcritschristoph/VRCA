import json
from scipy.spatial import distance
import sys
from Bio import SeqIO
from marker_genes import meta_marker
import operator
import argparse
import os.path
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
	result_string = ''
	for w in sorted(distances, key=distances.get):
		if count <= 3:
			result_string += "(" + w + ", " + str(distances[w]) + "), "
			count += 1
		else:
			break
	return result_string

def tetrat_compare(tetramers1, results):
	#Normalize
	distances = {}
	for tets in results[0]:
		tetramers2 = results[0][tets]
		name = results[1][tets] + " (" + tets + ")"

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
		distances[name] = round(distance.euclidean(query_dat, subject_dat),5)
	return sorted(distances.items(), key=operator.itemgetter(1))

if __name__ == "__main__":

	__author__ = "Alex Crits-Christoph"
	parser = argparse.ArgumentParser(description='Predicts host(s) for a contig through tetranucleotide similarity..')
	parser.add_argument('-i','--input', help='Viral contig(s) fasta file',required=True)
	parser.add_argument('-a','--assembly',help='Metagenome assembly FASTA file', required=False)
	parser.add_argument('-p', '--program', help="Program to run: 'd' for matching a database, 'm' for matching metagenomic assembled contigs, and 'md' for matching reference data for species found in the metagenomic assembly", required=False)

	args = parser.parse_args()
	if args.input:
		contigs = args.input
		print "Contigs"
	else:
		print "ERROR: No viral contig FASTA file was provided."
		sys.exit()


	if args.program:
		program = args.program
		if program == 'm' or program == 'md':
			if args.assembly:
				assembly = args.assembly
			else:
				print "ERROR: for programs 'm' and 'md' you must include a FASTA metagenome assembly using -a or --assembly."
				sys.exit()
		elif program == 'd':
			if args.assembly:
				print "Warning: ignoring metagenome assembly and running database-matching only."
	else:
		print "No program specified. Running database-matching only."

	if args.assembly:
		assembly = args.assembly
	else:
		print "ERROR: No viral contig FASTA file was provided."
		sys.exit()    
	if program == 'md' or program =='d':
		if os.path.isfile('./tetramer_database.dat'):
			database = 'tetramer_database.dat'
		else:
			print "ERROR: could not find tetramer database, and database-matching program was selected. Check for tetramer_database.dat in the local directory."
			sys.exit()
	else:
		database = 'tetramer_database.dat'
	#Calculate tetranucleotide frequencies for query
	handle = open(contigs, "rU")
	records = list(SeqIO.parse(handle, "fasta"))
	handle.close()

	metagenome = assembly
	if program == 'md' or program =='d':
		print "Loading database..."
		with open(database) as data_file:
			data = json.load(data_file)

	if program == 'm' or program == 'md':
		print "Finding marker contigs..."
		results = meta_marker.find_markers(metagenome)
	
	if program == 'md' or program == 'd':
		for seqs in records:
			print "Calculating tetranucleotide frequencies for query..."
			tetramers = calc_tetra(seqs)

			#Compare with database
			
			print "Comparing query tetranucleotide frequencies to database..."
			result_string = read_data(tetramers, data)
			print "*************** RESULTS ***************"
			print "Viral contig: Top database matches"
			print seqs.id + ": " + result_string 

	if program == 'm':
		for seqs in records:
			#Calculate Tetranucleotide frequencies for query
			print ""
			print "Calculating tetranucleotide frequencies for query..."
			tetramers = calc_tetra(seqs)

			#Compare with contigs
			print "Comparing query tetranucleotide frequencies to all marker contigs..."
			print "*************** RESULTS ***************"
			print "Viral contig: Top matches [blast_match (contig name)]" 
			print seqs.id + ": " + str(tetrat_compare(tetramers, results)[:3]).replace("]","").replace("[","")