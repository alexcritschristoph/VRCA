'''
Alex Crits-Christoph
License: GPL3
Finds Circular Contigs in Metagenomes using 2 different methods.
'''
import sys, getopt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import tempfile
import subprocess
import fractions

#Build dictionary of all contigs in input
def import_contigs(file_name):
    print "Importing contigs"
    handle = open(file_name, "rU")
    records = list(SeqIO.parse(handle, "fasta"))
    return records

#Run lastz on each contig, check for self-alignment at ends.
def self_align(seqs, read_length):
    print "Self-aligning contigs using lastz"
    print "Self-aligned contigs: "
    #Create temporary sequence file for lastz
    positive_self_aligns = []
    for seq in seqs:
        #Take last cutoff
        seq_part1 = str(seq.seq[0:read_length])
        seq_part2 = str(seq.seq[len(seq.seq)-read_length:len(seq.seq)])
        combined = seq_part1 + seq_part2
        f = tempfile.NamedTemporaryFile(delete=False)
        f.write(">" + seq.id + "\n")
        f.write(str(combined))
        f.seek(0)
        f.close()

        #run lastz
        output = subprocess.check_output(["lastz", f.name, "--self", "--notrivial", "--nomirror", "--format=general-:start1,end1,start2,end2,score,strand1,strand2,identity,length1"]);
        results = output.split("\n")
        for result in results:
            if result != '':
                start1 = result.split()[0]
                start2 = result.split()[2]
                end2 = result.split()[3]
                strand1 = result.split()[5]
                strand2 = result.split()[6]
                identity = result.split()[7]
                length = int(result.split()[9])
                if strand1 == strand2 and length > 0.5 * read_length and float(fractions.Fraction(identity)) > 0.9:
                    if int(start1) < 5 and int(end2) > read_length*2 * 0.9: 
                        print seq.id
			print result
			if seq.id not in positive_self_aligns:
                            positive_self_aligns.append(seq.id)
    return positive_self_aligns

if __name__ == "__main__":

    #Start with argument parsing
    try:
        opts, args = getopt.getopt(sys.argv[1:],"i:l:")
    except getopt.GetoptError:
        print 'find_circular.py -i input_assembly.fasta -l 100 '
        sys.exit(2)

    assembly = ''
    read_length = 0
    insert_size = 0
    output_dir = ''
    for opt, arg in opts:
        if opt == "-i":
            assembly = arg
        if opt == "-l":
            read_length = arg

    if assembly == '':
        print "ERROR: No contig FASTA file was provided."
        sys.exit()
    if read_length == 0:
        print "No read length specified. Using 100 bp."
    	read_length = 100

    #import contigs
    contigs = import_contigs(assembly)
    contigs_copy = contigs


    #check contigs for self-alignment overlap using lastz 
    self_aligned = self_align(contigs, read_length)
    
    #Write to file
    print "Writing to files..."
    seqs = {}
    for seq_r in contigs_copy:
        seqs[seq_r.id] = seq_r.seq

    f = open(assembly + '_circular_lastz.fna', 'a+')
    for contig in self_aligned:
	f.write(">" + str(contig) + "\n")
	f.write(str(seqs[contig]) + "\n")
    f.close()

    #Clean up
    print "Found " + str(len(self_aligned)) + " putative circular contigs with lastz."
    print "Completed. Output stored in " + assembly + '_circular_lastz.fna.'
