import sys
import numpy
from collections import defaultdict
contigs = []
genes = {}
introns = {}
confs = {}
strands = {}

with open(sys.argv[1], "r") as f:
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
			#print str(gene_size) + "," +  str(intron_length) + "," +  str(conf) + "," +  str(strandedness)

#print introns
for contig in contigs:
	if len(genes[contig]) > 20:
		print contig +"," + str(float(sum(genes[contig])) / len(genes[contig])) +  "," + str(float(sum(introns[contig])) / len(introns[contig])) + "," + str(float(sum(confs[contig])) / len(confs[contig])) + "," + str(float(sum(strands[contig])) / len(strands[contig]))
#print confs
