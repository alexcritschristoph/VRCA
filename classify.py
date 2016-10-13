import sys
import numpy
import subprocess
import os
import cPickle as pickle
from sklearn.ensemble import RandomForestClassifier
from collections import defaultdict
import os.path
from Bio import SeqIO
import argparse

def main(input_file, prodigal_path, min_contig):
	contigs = {}

	#Check that prodigal is installed
	try: 
	    subprocess.call(["which", prodigal_path])
	except: 
	    print "Error: Prodigal is not installed and available with the specified path."
	    sys.exit(1)

	#Check that the input file exists
	if not os.path.isfile(input_file):
		print "Error: The input file does not appear to exist."
		sys.exit(1)

	#Remove old output
	with open(os.devnull, 'wb') as devnull:
	    subprocess.check_call(['rm','-rf','./' + input_file.split("/")[-1] + '_output/'], stdout=devnull, stderr=subprocess.STDOUT)

	#Make new directory
	with open(os.devnull, 'wb') as devnull:
	    subprocess.check_call(['mkdir','./' + input_file.split("/")[-1] + '_output/'], stdout=devnull, stderr=subprocess.STDOUT)

	#Keep only contigs > min_length
	print "[1] Excluding contigs less than " + str(min_contig) + " bp in length."
	handle = open(input_file, "rU")
	l = SeqIO.parse(handle, "fasta")

	new_input = './' + input_file.split("/")[-1] + '_output/' + input_file.split("/")[-1] + "_filtered.fa"
	f = open(new_input, 'a+')
	contig_list = {}
	for s in l:
		if len(s.seq) >= int(min_contig):
			f.write(">" + s.id + "\n")
			f.write(str(s.seq) + "\n")
			contig_list[s.id] = s.seq
	f.close()

	#Run prodigal
	print "[2] Predicting genes and proteins with Prodigal (this could take a while)..."
	sts = subprocess.Popen(prodigal_path + ' -p meta -f gff -o ./'  + input_file.split("/")[-1] + '_output/temp.gff -i ' + new_input + ' -q', shell=True).wait()


	print "[3] Generating features for machine learning..."
	contigs = []
	genes = {}
	introns = {}
	confs = {}
	strands = {}

	with open('./' + input_file.split("/")[-1] + '_output/temp.gff', "r") as f:
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

	f = open('./' + input_file.split("/")[-1] + '_output/temp.dat','w+')
	for contig in contigs:
		if len(genes[contig]) > 10:
			f.write(str(float(sum(genes[contig])) / len(genes[contig])) +  "," + str(float(sum(introns[contig])) / len(introns[contig])) + "," + str(float(sum(confs[contig])) / len(confs[contig])) + "," + str(float(sum(strands[contig])) / len(strands[contig])) + "\n")


	f.close()

	print "[4] Loading classifier into memory..."
	try:
		forest_fit = pickle.load( open( "./data/classifier.p", "rb" ) )
	except:
		print "Error: could not load the classifier. Make sure it can be found in ./data/classifier.p."
		sys.exit(1)

	print "[5] Loading feature data..."
	testdata = numpy.loadtxt('./' + input_file.split("/")[-1] + '_output/temp.dat', delimiter=",")

	print "[6] Classifing features..."
	classified = list(forest_fit.predict(testdata))
	j = 0
	print "Viral contigs:"
	viral_contigs = []
	results_string = ''
	for i in classified:
		if i == 1:
			results_string +=  str(contigs[j]) + ","
			viral_contigs.append(contigs[j])
		j += 1

	results_string = results_string.rstrip(",")
	print results_string 
	print "[7] Saving results to ./" + input_file.split("/")[-1] + '_output/viral.fna'
	f = open('./' + input_file.split("/")[-1] + '_output/viral.fna', 'a+')
	for contig in viral_contigs:
		f.write(">" + str(contig_list[contig]) + "\n")
		f.write(str(contig_list[contig]))

	print "[8] Complete. Cleaning up and exiting..."

if __name__ == '__main__':

    __author__ = "Alex Crits-Christoph"

    parser = argparse.ArgumentParser(description='Classifies viral metagenomic contigs using a random forest classifier built on public datasets. Features are gene length, intergenic space length, strandedness, and prodigal calling gene confidence.')
    parser.add_argument('-i','--input', help='Input assembly filename',required=True)
    parser.add_argument('-m', '--min_contig_size', help='Minimum contig size to use (default: 10 kbp)', required=False)
    parser.add_argument('-p', '--prodigal_path', help='Path to prodigal (default: prodigal)', required=False)

    args = parser.parse_args()
    if args.min_contig_size == None:
    	args.min_contig_size = 10000
    if args.prodigal_path == None:
    	args.prodigal_path = 'prodigal'
	
	main(args.input, args.prodigal_path, args.min_contig_size)