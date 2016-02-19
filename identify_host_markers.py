import json
from scipy.spatial import distance
import sys
from Bio import SeqIO
from marker_genes import meta_marker
import operator
import argparse
import os.path
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

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
	real_names = {}
	for tets in results[0]:
		tetramers2 = results[0][tets]
		name = results[1][tets] + " (" + tets + ")"
		real_names[results[1][tets] + " (" + tets + ")"] = tets

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
	return [real_names, sorted(distances.items(), key=operator.itemgetter(1))]

def visualize(data):
	print "visualizing results..."
	tetramer_array = []
	names = data[0][1]
	sizes = data[0][2]
	contig_count = 0
	name_pos = []
	for t in sorted(data[0][0].keys()):
		name_pos.append(t)
		contig_count += 1
		temp = []
		for tet in sorted(data[0][0][t].keys()):
			temp.append(data[0][0][t][tet])
		tetramer_array.append(temp)

	for t in sorted(data[1].keys()):
		temp = []
		for tet in sorted(data[1][t].keys()):
			temp.append(data[1][t][tet])
		tetramer_array.append(temp)

	tetramers_np = np.array(tetramer_array)
	pca = PCA(n_components=2)
	fit = pca.fit(tetramers_np).transform(tetramers_np)

	#Contig sizing
	for contig in sizes:
		sizes[contig] = float(sizes[contig]) / float(max(sizes.values())) * 250 + 35
	sizes_list = []
	for contig in sorted(sizes.keys()):
		sizes_list.append(sizes[contig])

	#Contig colors
	color_assigned = {}
	color_list = []
	color_brewer = ['#1f78b4', '#33a02c', '#e31a1c', '#ff7f00', '#6a3d9a', '#b15928', '#a6cee3', '#b2df8a', '#fb9a99', '#fdbf6f', '#cab2d6', '#ffff99']
	i = 0
	#assign colors
	most_common_names = Counter(names.values()).most_common()
	for name in most_common_names:
		if len(name) > 1:
			n = name[0]
			if n not in color_assigned:
				if i < len(color_brewer)-1:
					color_assigned[n] = color_brewer[i]
					i += 1
				else:
					color_assigned[n] = color_brewer[i]
		else:
			color_assigned[n] = '#ffff99'
	#create color list
	for name in sorted(names.keys()):
		color_list.append(color_assigned[names[name]])

	#plot it
	plt.scatter(fit[0:contig_count,0], fit[0:contig_count:,1], s=sizes_list,  c=color_list, marker= 'o', alpha=0.9)
	plt.scatter(fit[contig_count:,0], fit[contig_count:,1], marker= 'x', s=50, c="#ef1a1a", alpha=1)

	#Add labels
		#Add labels
	positions = {}
	i = 0 
	avg_size = sum(sizes.values()) / len(sizes.values())
	for name in sorted(names.keys()):
		if len(names.keys()) > 5 and sizes[name] > avg_size:
			positions[names[name]] = [fit[i,0], fit[i,1]]
		i += 1

	#Plot virus - closest host lines
	i = 0
	for t in sorted(data[1].keys()):
		#get pos of name in fit
		name = name_pos.index(data[2][t])
		#plot line
		plt.plot([fit[contig_count+i,0], fit[name,0]], [fit[contig_count+i,1], fit[name,1]], alpha=0.5)

		plt.annotate(
        t, 
        xy = (fit[contig_count+i,0], fit[contig_count+i,1]), xytext = (0, -20),
        textcoords = 'offset points', ha = 'right', va = 'bottom',
        fontsize = 8,
        bbox = dict(boxstyle = 'round,pad=0.2', fc = 'white', alpha = 0.5),
        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

		i += 1

	for name in positions.keys():
		x = positions[name][0]
		y = positions[name][1]
		
		plt.annotate(
        name, 
        xy = (x, y), xytext = (-30, 30),
        textcoords = 'offset points', ha = 'right', va = 'bottom',
        fontsize = 10,
        bbox = dict(boxstyle = 'round,pad=0.2', fc = 'white', alpha = 0.5),
        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))


	plt.show()

if __name__ == "__main__":

	__author__ = "Alex Crits-Christoph"
	parser = argparse.ArgumentParser(description='Predicts host(s) for a contig through tetranucleotide similarity..')
	parser.add_argument('-i','--input', help='Viral contig(s) fasta file',required=True)
	parser.add_argument('-a','--assembly',help='Metagenome assembly FASTA file', required=False)
	parser.add_argument('-p', '--program', help="Program to run: 'd' for matching a database, 'm' for matching metagenomic assembled contigs, and 'md' for matching reference data for species found in the metagenomic assembly", required=False)
	parser.add_argument('-v', '--visualize', help="Visualizes a PCA of tetranucleotide frequencies of host contigs and viral contigs.", required=False)

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
		#Compare with database
		print "Comparing query tetranucleotide frequencies to database..."
		print "*************** RESULTS ***************"
		print "Viral contig: Top database matches"
		for seqs in records:
			tetramers = calc_tetra(seqs)
			result_string = read_data(tetramers, data)
			print seqs.id + ": " + result_string 

	if program == 'm' or program == 'md':
		#Compare with contigs
		print "Comparing query tetranucleotide frequencies to all marker contigs..."
		print "*************** RESULTS ***************"
		print "Viral contig: Top matches [blast_match (contig name)]" 
		combined_results = [results, {}, {}]
		for seqs in records:
			tetramers = calc_tetra(seqs)
			tetrats = tetrat_compare(tetramers, results)
			top_three = tetrats[1][:3]
			print seqs.id + ": " + str(top_three).replace("]","").replace("[","")
			combined_results[1][seqs.id] = tetramers
			combined_results[2][seqs.id] = tetrats[0][top_three[0][0]]
		if args.visualize:
			visualize(combined_results)