#!/usr/bin/env python
'''
Alex Crits-Christoph
License: GPL3
Finds Circular Contigs in assembled metagenomes by checking for start/end forward alignment.
'''
import sys
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import tempfile
import subprocess
import fractions

#Build dictionary of all contigs in input
def import_contigs(file_name, min_contig_size):
    print "Importing contigs"
    handle = open(file_name, "rU")
    records = list(SeqIO.parse(handle, "fasta"))
    records_filtered = []
    for seq in records:
        if len(seq.seq) >= min_contig_size:
            records_filtered.append(seq)
    return records_filtered

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
        f = tempfile.NamedTemporaryFile(delete=True)
        f.write(">" + seq.id + "\n")
        f.write(str(combined))
        f.seek(0)
        
        #run lastz
        output = subprocess.check_output(["lastz", f.name, "--self", "--notrivial", "--nomirror", "--format=general-:start1,end1,start2,end2,score,strand1,strand2,identity,length1"]);
        f.close()
        results = output.split("\n")
        for result in results:
            if result != '':
                start1 = result.split()[0]
                end1 = result.split()[1]
                start2 = result.split()[2]
                end2 = result.split()[3]
                strand1 = result.split()[5]
                strand2 = result.split()[6]
                identity = result.split()[7]
                length = int(result.split()[9])
                if strand1 == strand2 and length > 0.4 * read_length and float(fractions.Fraction(identity)) > 0.95:
                    if int(start1) < 5 and int(start2) > read_length and int(end1) < read_length and int(end2) > read_length*2 * 0.9: 
                        print seq.id
                        print result
                        if seq.id not in positive_self_aligns:
                            positive_self_aligns.append(seq.id)
    return positive_self_aligns

if __name__ == "__main__":

    __author__ = "Alex Crits-Christoph"

    parser = argparse.ArgumentParser(description='Finds circular contigs using lastz.')
    parser.add_argument('-i','--input', help='Input assembly filename',required=True)
    parser.add_argument('-l','--read_length',help='Read length (default: 101 bp) [this is what will be checked for overlap]', required=False)
    parser.add_argument('-m', '--min_contig_size', help='Minimum contig size to check for (default: 3 kbp)', required=False)

    args = parser.parse_args()


    assembly = ''
    read_length = 101
    min_contig = 3000
    output_dir = ''

    if args.min_contig_size:
        min_contig = int(args.min_contig_size)

    if args.input:
        assembly = args.input
    else:
        print "ERROR: No contig FASTA file was provided."
        sys.exit()


    if args.read_length:
        read_length = int(args.read_length)
           
    #Remove old temporary files if they exist.
    print "Removing old temp files [Nope, I'm not -Matt]"
    #os.system('rm ' + assembly + '.* > /dev/null 2>&1')
    #os.system('rm ' + assembly + '_* > /dev/null 2>&1')

    #import contigs [all this is doing is filtering by size]
    contigs = import_contigs(assembly, min_contig)
    contigs_copy = contigs

    self_aligned = []
    #check contigs for self-alignment overlap using lastz 
    self_aligned = self_align(contigs, read_length)
    
    paired_results = []

    #Write to file
    print "Writing to files..."
    seqs = {}
    if len(self_aligned) > 0:
        for seq_r in contigs_copy:
            seqs[seq_r.id] = seq_r.seq
        f = open(assembly.split("/")[-1] + '_circular.fna', 'a+')
        for contig in self_aligned:
            f.write(">" + str(contig) + "\n")
            f.write(str(seqs[contig]) + "\n")
        f.close()
    else:
        print "Error: No circular contigs were found."
        sys.exit(1)
    #Clean up
    #os.system('rm ' + assembly + '.* > /dev/null 2>&1')
    #os.system('rm ' + assembly + '_output* > /dev/null 2>&1')
    #os.system('rm ' + assembly + '_overlap* > /dev/null 2>&1')
    print "Found " + str(len(self_aligned)) + " putative circular contigs."
    print "Completed. Output stored in " + assembly.split("/")[-1] + '_circular.fna.'

