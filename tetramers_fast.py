from Bio import SeqIO
import sys
import glob
import khmer

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
	print tetramers
if __name__ == "__main__":

	handle = open(sys.argv[1], "rU")
	records = list(SeqIO.parse(handle, "fasta"))
	handle.close()
	calc_tetra(records)