'''
Alex Crits-Christoph
License: GPL3
Finds Circular Contigs in Metagenomes using 2 different methods.
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

#Analyze bowtie alignment output
def calculate_paired_circle(samfile, read_length, insert_size):
    print "Finding circular tips..."
    import pysam
    samfile = pysam.AlignmentFile(samfile,'rb')
    circular = 0
    not_c = 0
    true_circle = []
    ids = []
    all_reads = []
    
    for r in samfile.references:
        circle1 = False
        circle_reads = 0
        total = 0
        for read in samfile.fetch(r):
            if read.query_alignment_length > 0.9*read_length and read.reference_end < read_length+insert_size and read.next_reference_start > (read_length+insert_size) +5:
                    circle1 = True
                    circle_reads += 1
            total += 1

        if circle1 == True:
            circular += 1
            true_circle.append(circle_reads)
            all_reads.append(total)
            ids.append(r)
        else:
                not_c += 1

    i = 0 
    final_circle = []
    print "Identified Circular Elements: "
    for read in all_reads:
        percent = float(true_circle[i]) / read      
        if percent > 0.01:
                print ids[i] + ":" + str(round(percent,3))
                final_circle.append(ids[i])
 	else:
		print "Didn't make cutoff of 1% of reads bridging: " + ids[i] + ":" + str(round(percent,3))
        i += 1
    print "Stats: "
    print "putative circular elements: " + str(circular)
    print "no evidence for circularity: " + str(not_c)
    return final_circle

def create_edges(contigs, assembly, overlap_size):
    print "Extracting contig edges..."
    overlap_seqs = {}

    for record in contigs:
        overlap_seqs[record.id] = record.seq[-overlap_size:] + record.seq[:overlap_size]

    #Writing to file...
    f = open(assembly + ".temp", 'a+')

    for seq in overlap_seqs.keys():
       f.write(">" + str(seq) + "\n")
       f.write(str(overlap_seqs[seq]) + "\n")
    f.close()

    return assembly + '.temp'

def run_bowtie(assembly, reads_1, reads_2):
    print "Creating bowtie database..."
    #Run bowtie-build on overlap seqs
    os.system('bowtie2-build ' + assembly + '.temp ' + assembly + "_overlap -q")
 
    print "Mapping reads to overlaps with bowtie... (this could take up to an hour)"

    #Run bowtie on paired reads on overlap seqs.
    print 'bowtie2 -x ' + assembly + '_overlap -1 ' + reads_1 + ' -2 ' + reads_2 + ' -S ' + assembly + '_temp.sam --no-unal --no-mixed --no-discordant -p 8 -I ' + str(insert_size / 2) + ' -X ' + str(insert_size * 2)
    os.system('bowtie2 -x ' + assembly + '_overlap -1 ' + reads_1 + ' -2 ' + reads_2 + ' -S ' + assembly + '_temp.sam --no-unal --no-mixed --no-discordant -p 8 -I ' + str(insert_size / 2) + ' -X ' + str(insert_size * 2))

    #Find those with overlapping matches
    print "Converting sam to bam..."
    os.system('samtools view -bS ' + assembly + '_temp.sam > ' + assembly + '_output.bam')
    print "Sorting and indexing bam..."
    os.system('samtools sort ' + assembly + '_output.bam ' + assembly + '_output_sorted.bam')
    os.system('samtools index ' + assembly + '_output_sorted.bam.bam')
    return assembly + '_output_sorted.bam.bam'

if __name__ == "__main__":

    __author__ = "Alex Crits-Christoph"
    print sys.argv

    parser = argparse.ArgumentParser(description='This is a demo script by nixCraft.')
    parser.add_argument('-i','--input', help='Input assembly filename',required=True)
    parser.add_argument('-r1','--reads1',help='Paired reads file 1', required=False)
    parser.add_argument('-r2','--reads2',help='Paired reads file 2', required=False)
    parser.add_argument('-l','--read_length',help='Read length (default: 101 bp)', required=False)
    parser.add_argument('-s','--insert_size',help='Insert Size (default: 150 bp)', required=False)
    parser.add_argument('-m', '--min_contig_size', help='Minimum contig size to check for (default: 3 kbp)', required=False)
    parser.add_argument('-p', '--program',help="Program - 'lastz', 'bowtie', or 'both'(default)", required=False)

    args = parser.parse_args()


    assembly = ''
    reads1 = ''
    reads2 = ''
    read_length = 101
    insert_size = 150
    min_contig = 3000
    output_dir = ''

    if args.min_contig_size:
        min_contig = int(args.min_contig_size)

    if args.input:
        assembly = args.input
    else:
        print "ERROR: No contig FASTA file was provided."
        sys.exit()

    if args.program and args.program != "lastz" and args.program != "bowtie" and args.program != "both":
        print "ERROR: invalid program type. Use -p lastz, -p bowtie, or -p both."
        sys.exit()

    if not args.program or (args.program and args.program != 'lastz'):
        if args.reads1:
            reads1 = args.reads1
        else:
            print "ERROR: Paired Reads file #1 not found. To run without reads file, use -p lastz."
            sys.exit()

        if args.reads2:
            reads2 = args.reads2
        else:
            print "ERROR: Paired Reads file #2 not found. To run without reads file, use -p lastz."
            sys.exit()
    if args.program == 'lastz':
        program = 'lastz'
    elif args.program == 'bowtie':
        program = 'bowtie'
    else:
        program = 'both'

    if args.read_length:
        read_length = int(args.read_length)
    if args.insert_size:
        insert_size = int(args.insert_size)
           
    #Remove old temporary files if they exist.
    print "Removing old temp files"
    os.system('rm ' + assembly + '.* > /dev/null 2>&1')
    os.system('rm ' + assembly + '_* > /dev/null 2>&1')

    #import contigs
    contigs = import_contigs(assembly, min_contig)
    contigs_copy = contigs

    overlap_size = read_length + insert_size
    self_aligned = []
    #check contigs for self-alignment overlap using lastz 
    if program != 'bowtie':
        self_aligned = self_align(contigs, read_length)
    
    #Write contig edges to file
    contig_edges = create_edges(contigs, assembly, overlap_size)

    paired_results = []
    #Run bowtie
    if program != 'lastz':
        sam_file = run_bowtie(assembly, reads1, reads2)
        #Analyze results
        paired_results = calculate_paired_circle(sam_file, read_length,insert_size)

    #Write to file
    print "Writing to files..."
    seqs = {}
    for seq_r in contigs_copy:
        seqs[seq_r.id] = seq_r.seq
    if program != 'bowtie':
        f = open(assembly + '_circular_lastz.fna', 'a+')
        for contig in self_aligned:
     	    f.write(">" + str(contig) + "\n")
	    f.write(str(seqs[contig]) + "\n")
        f.close()
    if program != 'lastz':
        f = open(assembly + '_circular_bowtie.fna', 'a+')
        for contig in paired_results:
	    f.write(">" + str(contig) + "\n")
	    f.write(str(seqs[contig]) + "\n")
        f.close()

    if program == 'both':
        consensus = list(set(paired_results).intersection(self_aligned))

        f = open(assembly + '_circular_consensus.fna', 'a+')
        for contig in consensus:
            f.write(">" + str(contig) + "\n")
            f.write(str(seqs[contig]) + "\n")
        f.close()
    
    #Clean up
    os.system('rm ' + assembly + '.* > /dev/null 2>&1')
    os.system('rm ' + assembly + '_output* > /dev/null 2>&1')
    os.system('rm ' + assembly + '_overlap* > /dev/null 2>&1')
    print "Found " + str(len(self_aligned)) + " putative circular contigs with lastz and " + str(len(paired_results)) + " putative circular contigs using bowtie2. " + str(len(consensus)) + " contigs were found by both methods."
    print "Completed. Output stored in " + assembly + '_circular_lastz.fna, ' + assembly + '_circular_bowtie.fna, and ' + assembly + "_circular_consensus.fna."

