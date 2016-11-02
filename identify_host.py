'''
Alex Crits-Christoph
License: GPL3
Identifies closest GenBank prokaryote genomes by euclidean distance between tetranucleotide frequencies of query sequences.

'''

import json
from scipy.spatial import distance
import sys
from Bio import SeqIO
from marker_genes import meta_marker
import operator
import argparse
import os.path
from sklearn import manifold
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter


##Calculates the tetramer counts for an input sequence
def calc_tetra(seqs):
	#Create tetramers dict
	tetramers = {}
	for a in ['A', 'C', 'G', 'T']:
		for b in ['A', 'C', 'G', 'T']:
			for c in ['A', 'C', 'G', 'T']:
				for d in ['A', 'C', 'G', 'T']:
					tetramers[a+b+c+d] = 0

	#Count tetramers across sequence in a 4 bp sliding window
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

	#Return tetramers dictionary
	return tetramers

##
##
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

def visualize(data, file_name):
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
	seed = np.random.RandomState(seed=3)
	mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=seed,
                   dissimilarity="euclidean", n_jobs=1)
	fit = mds.fit_transform(tetramers_np)

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
	
	#Output graphical data to file
	f = open('./pca_data.txt', 'w')
	f.write("Host contigs:\n")
	f.write("Contig name, BLAST hit, x coord, y coord\n")
	i = 0
	for t in sorted(data[0][0].keys()):
		f.write(str(t) + "," + str(names[t]) + "," + str(fit[i,0]) + "," + str(fit[i,1]) + "\n")
		i += 1
	f.write("Viral contigs:\n")
	f.write("Contig name, closest host contig, closest host name, x coord, y coord\n")

	for t in sorted(data[1].keys()):
		name = data[2][t]
		
		f.write(str(t) + "," + str(name) + "," + str(names[name]) + "," + str(fit[i,0]) + "," + str(fit[i,1]) + "\n")
		i += 1
	f.close()
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

       #		plt.annotate(
       # t, 
       # xy = (fit[contig_count+i,0], fit[contig_count+i,1]), xytext = (0, -20),
       # textcoords = 'offset points', ha = 'right', va = 'bottom',
       # fontsize = 8,
       # bbox = dict(boxstyle = 'round,pad=0.2', fc = 'white', alpha = 0.5),
       # arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

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

	plt.savefig(file_name)
	plt.show()


## Main Function
if __name__ == "__main__":

	## Handle arguments.
	__author__ = "Alex Crits-Christoph"
	parser = argparse.ArgumentParser(description='Predicts host(s) for a contig through tetranucleotide similarity.')
	parser.add_argument('-i','--input', help='Viral contig(s) fasta file',required=True)
	parser.add_argument('-a','--assembly',help='Cellular metagenome assembly FASTA file (to search for host marker genes).', required=False)
	parser.add_argument('-v', '--visualize', help="Visualizes NMDS of tetranucleotide frequencies for host contigs and viral contigs.", required=False, action='store_true')
	parser.add_argument('-b', '--blast_path', help="Path to the BLASTP executable (default: blastp).", required=False)
	parser.add_argument('-hm', '--hmmsearch_path', help="Path to the hmmsearch executable (default: hmmsearch).", required=False)
	parser.add_argument('-o', '--output', help="Path to the directory to store output (will create if does not exist).", required=False)

	args = parser.parse_args()

	#Check to see if output directory exists
	if args.output:
		if os.path.isdir(args.output):
			output_dir = args.output
		else:
			os.system("mkdir " + args.output)
			if os.path.isdir(args.output):
				output_dir = args.output
			else:
				print "[ERROR] Could not create output directory."
				sys.exit(1)

	else:
		output_dir = './'

	blast_path = 'blastp'
	if args.blast_path:
		blast_path = args.blast_path

	hmmsearch_path = 'hmmsearch'
	if args.hmmsearch_path:
		hmmsearch_path = args.hmmsearch_path

	if args.input:
		if os.path.isfile(args.input):
			contigs = args.input
		else:
			print "[ERROR] The provided viral contig FASTA file was not found"
			sys.exit()
	else:
		print "[ERROR] No viral contig FASTA file was provided"
		sys.exit()

	database = './data/tetramer_database.dat'
	if os.path.isfile('./data/tetramer_database.dat'):
		print "[SETUP] Found tetramer database file, using " + database + "."
	else:
		print "[ERROR] could not find tetramer database. Check for ./data/tetramer_database.dat in the local directory"
		sys.exit()

	if args.assembly:
		if os.path.isfile(args.input):
			print "[SETUP] Assembly included. Will search for host marker genes in the assembly and compare viral contigs to genomes of identified hosts"
			metagenome = args.assembly
		else:
			print "[ERROR] The provided metagenome assembly contig FASTA file was not found"
			sys.exit()
	else:
		print "[SETUP] No assembly file included. Will compare viral contigs to all host genomes in GenBank"

	#Read in viral queries
	handle = open(contigs, "rU")
	records = list(SeqIO.parse(handle, "fasta"))
	handle.close()
	
	print "[EXEC] Loading database"
	with open(database) as data_file:
		data = json.load(data_file)

	if args.assembly:
		print "[EXEC] Searching for marker contigs in metagenome assembly"
		results = meta_marker.find_markers(metagenome, blast_path, hmmsearch_path, output_dir)
	
	# #Compare with database
	# print "[EXEC] Comparing query tetranucleotide frequencies to entire database"
	# f = open(output_dir.rstrip("/") + "/" + contigs.split("/")[-1] + "_genbank.txt", 'w+')
	# f.write("Viral contig name: (Closest GenBank match, Euclidean Distance)\n")
	# for seqs in records:
	# 	tetramers = calc_tetra(seqs)
	# 	result_string = read_data(tetramers, data)
	# 	f.write(seqs.id + ": " + result_string + "\n")
	# f.close()
	
	if args.assembly:
		#Compare with contigs
		print "[EXEC] Comparing query tetranucleotide frequencies to contigs with marker proteins"
		f = open(output_dir.rstrip("/") + "/" + contigs.split("/")[-1].split(".")[0] + "_markercontigs.txt", 'w+')
		f.write("Viral contig name: (Closest contig identity (contig name), Euclidean Distance to marker contig)\n")
		combined_results = [results, {}, {}]
		for seqs in records:
			tetramers = calc_tetra(seqs)
			tetrats = tetrat_compare(tetramers, results)
			top_three = tetrats[1][:3]
			f.write(seqs.id + ": " + str(top_three).replace("]","").replace("[","") + "\n")
			combined_results[1][seqs.id] = tetramers
			combined_results[2][seqs.id] = tetrats[0][top_three[0][0]]
		f.close()

		#Compare with closest BLAST matches in the database
		print "[EXEC] Comparing query tetranucleotide frequencies to genomes of species identified in metagenome"
		f = open(output_dir.rstrip("/") + "/" + contigs.split("/")[-1].split(".")[0] + "_markergenbank.txt", 'w+')
		f.write("Viral contig name: (Identified species (marker contig name), Euclidean Distance to GenBank genome)\n")
		
		#Create subset of genomes which were found in this assembly
		genomes_in_metagenome = {}
		species_to_contigs = {}
		for contig in results[1]:
			genus = results[1][contig]
			for species in data.keys():
				if genus in species:
					genomes_in_metagenome[species + " (" + contig + ")"] = data[species]
					species_to_contigs[species] = contig
		for seqs in records:
			tetramers = calc_tetra(seqs)
			result_string = read_data(tetramers, genomes_in_metagenome)
			f.write(seqs.id + ": " + result_string + "\n")
		f.close()

		#Visualization
		print "[EXEC] Visualizing results"
		if args.visualize:
			visualize(combined_results, output_dir.rstrip("/") + "/" + contigs.split("/")[-1].split(".")[0] + ".png")
