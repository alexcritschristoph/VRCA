'''
Finds CRISPR spacers in a metagenome; then looks to see if those spacers have good BLAST matches to other contigs in that metagenome.
The implication is the CRISPR spacer will be on a host contig and the BLAST match on the other contig could be a phage it is targeting.
'''

import sys
import argparse
import subprocess
import os
from Bio import SeqIO
from Bio.SeqUtils import GC

#convenient code that checks path for minced and blastn.
def which(program):
	import os
	def is_exe(fpath):
		return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

	fpath, fname = os.path.split(program)
	if fpath:
		if is_exe(program):
			return program
	else:
		for path in os.environ["PATH"].split(os.pathsep):
			path = path.strip('"')
			exe_file = os.path.join(path, program)
			if is_exe(exe_file):
				return exe_file

	return None

##Check for valid input 

if __name__ == "__main__":
	__author__ = "Alex Crits-Christoph"
	parser = argparse.ArgumentParser(description='Identifies CRISPR arrays in an assembled metagenome, and matches identified spacers to other contigs.')
	parser.add_argument('-i','--input', help='Input assembly filename (FASTA format)',required=True)

	args = parser.parse_args()

	if args.input:
		assembly = args.input
	else:
		print("ERROR: No assembly file indicated. Use like: crispr_matches -i assembly.fna")
		sys.exit()

	if not which("minced"):
		print("ERROR: Please install the program minced in your path.")
		sys.exit()

	if not which("makeblastdb"):
		print("ERROR: Please install BLAST+ (makeblastdb) in your path.")
		sys.exit()		

	if not which("blastn"):
		print("ERROR: Please install BLAST+ (blastn) in your path.")
		sys.exit()		

	print("All requirements met. Creating new output directory...")
	os.system("rm -rf " + assembly + "_crispr_matches")
	os.system("mkdir " + assembly + "_crispr_matches")
	os.system("cp " + assembly + " ./" + assembly+ "_crispr_matches")

	filename = "./" + assembly+ "_crispr_matches/" + assembly
	print("Finding CRISPRs with minced...")
	output = subprocess.check_output(["minced", "-gff", "-spacers", filename], encoding='utf8')
#	lines = output.split("\n")


	#Read fasta file
	handle = open(filename, "rU")
	records_dict = {}
	records = SeqIO.parse(handle, "fasta")
	for record in records:
		records_dict[record.id] = str(record.seq)
	handle.close()

	crispr_contigs = []
	for line in output.split("\n"):
		if line.strip() != '' and not line.startswith("#"):
			contig = line.split()[0]
			crispr_contigs.append(contig)
	crisprs = open("./" + assembly+ "_crispr_matches/" + assembly + "_crisprs.fna", 'a+')
	for contig in crispr_contigs:
		crisprs.write(">" + contig + "\n")
		crisprs.write(records_dict[contig] + "\n")
	print(str(len(crispr_contigs)) + " contigs with CRISPRs found, stored in " + assembly+ "_crispr_matches/" + assembly + "_crisprs.fna")

	print("Creating blastdb for non-CRISPR contigs...")

	f = open("./" + assembly+ "_crispr_matches/" + assembly + '_nocrisprs.fna', 'a+')
	for contig in records_dict:
		if contig not in crispr_contigs:
			f.write(">" + contig + "\n")
			f.write(records_dict[contig] + "\n")
	f.close()

	output = subprocess.check_output("makeblastdb" + " -in " + "./" + assembly+ "_crispr_matches/" + assembly + '_nocrisprs.fna' + " -dbtype nucl", shell=True, encoding='utf8')

	print("BLASTing spacers against non-CRISPR contigs...")
	output = subprocess.check_output("blastn " + "-query " + "./" + assembly+ "_crispr_matches/" + assembly.split(".")[0] + "_spacers.fa" + " -db " + "./" + assembly+ "_crispr_matches/" + assembly + '_nocrisprs.fna' + " -outfmt 6", shell=True, encoding='utf8')
	
	print("Host/CRISPR contig\tHost contig length\tSpacer #\tViral contig\tViral contig length\tMatch PID\tMatch Length\tHost GC%\tViral GC%")
	for line in output.split("\n"):
		try:
			crispr_contig = line.split("_CRISPR")[0]
			crispr_spacer = line.split("spacer_")[1].split()[0]
			viral_contig = line.split("\t")[1]
			crispr_length = str(len(records_dict[crispr_contig]))
			viral_length = str(len(records_dict[viral_contig]))
			pid = line.split("\t")[2]
			length = line.split("\t")[3]
			host_gc = str(GC(records_dict[crispr_contig]))
			viral_gc = str(GC(records_dict[viral_contig]))
			print(crispr_contig + "\t" + crispr_length + "\t" + crispr_spacer + "\t" + viral_contig + "\t" + viral_length + "\t" + pid + "\t" + length + "\t" + host_gc + "\t" + viral_gc)
		except:
			pass
