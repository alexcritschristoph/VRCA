import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from random import randint

j = 1
import sys
for root, subdirs, files in os.walk(sys.argv[1]):
	for f in files:
		seq = os.path.join(root,f)
		seqr = SeqIO.read(open(seq), "fasta")
		limit=len(seqr.seq)
		#generate a random number of fragments from this genome
		#rand(1,5) for bacteria; rand(1,1) for viruses.
		#random size, greater than 10,000 bp and smaller than 0.75 of its genome)
		if limit > 10000:
			#Maximum is 100,000 bp or 0.75 * genome size
			max_s = int(0.75*limit)
			if max_s > 100000:
				max_s = 100000
			size = randint(1000,max_s)
			#random start- anywhere from 0 to genome_length - size
			start = randint(0,limit-size)
			end = start + size

			fragment = seqr.seq[start:end]

			print ">" + str(j)
			print fragment

		j += 1