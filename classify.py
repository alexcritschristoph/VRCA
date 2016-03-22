import sys
import numpy
import subprocess
import os

from collections import defaultdict
contigs = {}

with open(os.devnull, 'wb') as devnull:
    subprocess.check_call(['rm','-rf','./' + sys.argv[1] + '_output/'], stdout=devnull, stderr=subprocess.STDOUT)

with open(os.devnull, 'wb') as devnull:
    subprocess.check_call(['mkdir','./' + sys.argv[1] + '_output/'], stdout=devnull, stderr=subprocess.STDOUT)

print "[1] Predicting genes and proteins with Prodigal (this could take a while)..."
sts = subprocess.Popen('prodigal -a ./' + sys.argv[1] + '_output/temp.faa -p meta -f gff -o ./'  + sys.argv[1] + '_output/temp.gff -i ' + sys.argv[1] + ' -q', shell=True).wait()


print "[2] Generating features for machine learning..."
contigs = []
genes = {}
introns = {}
confs = {}
strands = {}

with open('./' + sys.argv[1] + '_output/temp.gff', "r") as f:
	previous_end = 0
	previous_strand = -1
	previous_contig = ""
	for line in f:
		if line.startswith("#"):
			previous_end = 0
			previous_strand = 0
		if not line.startswith("#"):
			#get gene size
			contig = line.split()[0]
			s = int(line.split()[3])
			e = int(line.split()[4])
			gene_size = int(e) - int(s)

			#get reading strand (-) or (+)
			stra = line.split()[6]
			if stra == "+":
				strand = 1.0
			else:
				strand = 0.0

			if contig != previous_contig:
				previous_contig = contig
				previous_end = e
				previous_strand = strand
				contigs.append(contig)
				genes[contig] = []
				introns[contig] = []
				confs[contig] = []
				strands[contig] = []
				continue

			#get intron length
			intron_length = s - previous_end
			if strand == previous_strand:
				strandedness = 1
			else:
				strandedness = 0

			previous_end = e
			previous_strand = strand

			#get conf
			conf = float(line.split("conf=")[1].split(";")[0])

			genes[contig].append(gene_size)
			introns[contig].append(intron_length)
			confs[contig].append(conf)
			strands[contig].append(strandedness)

f = open('./' + sys.argv[1] + '_output/temp.dat','w+')
for contig in contigs:
	if len(genes[contig]) > 10:
		f.write(str(float(sum(genes[contig])) / len(genes[contig])) +  "," + str(float(sum(introns[contig])) / len(introns[contig])) + "," + str(float(sum(confs[contig])) / len(confs[contig])) + "," + str(float(sum(strands[contig])) / len(strands[contig])) + "\n")


f.close()

print "[3] Loading classifier into memory..."
import cPickle as pickle
from sklearn.ensemble import RandomForestClassifier
forest_fit = pickle.load( open( "./data/classifier.p", "rb" ) )

print "[4] Loading feature data..."
testdata = numpy.loadtxt('./' + sys.argv[1] + '_output/temp.dat', delimiter=",")

print "[5] Classifing features..."
classified = list(forest_fit.predict(testdata))

print "[6] Saving results to ./output/viral.fna..."
j = 0
for i in classified:
	print str(contigs[j]) + ":" + str(i)
	j += 1

print "[7] Complete. Cleaning up and exiting..."
