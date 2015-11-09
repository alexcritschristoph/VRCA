from Bio import SeqIO
import sys
import glob
import khmer
import json

def calc_tetra(seqs):

	ktable = khmer.new_ktable(4)

	for seq in seqs:
		ktable.consume(str(seq.seq))

	tetramers = {}
	for a in ['A', 'C', 'G', 'T']:
		for b in ['A', 'C', 'G', 'T']:
			for c in ['A', 'C', 'G', 'T']:
				for d in ['A', 'C', 'G', 'T']:
					tetramers[a+b+c+d] = 0
	for t in tetramers.keys():
		tetramers[t] = int(ktable.get(t))
	return tetramers

if __name__ == "__main__":

	#For file in directory
	genbank_data = {}
	count = 0

	for filen in glob.glob('./genomes/*.fna'):
		print "selecting file"
		handle = open(filen, "rU")
		records = list(SeqIO.parse(handle, "fasta"))
		handle.close()
		name = filen.split("/")[-1].split(".fna")[0].replace("_"," ")

		print "calculating freqs..."

		genbank_data[name] = calc_tetra(records)
		
		count += 1
		print count

	f = open('./tetra.dat', 'a+')
	f.write(json.dumps(genbank_data))
	f.close
	print ("Success!")