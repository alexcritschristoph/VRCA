## VIral and Circular content from metAgenomes (VICA)

**`find_circular_fast.py`**

Quickly finds circular contigs in metagenome assemblies by looking for the overlap at the start and end of contigs. Please run `pip install edlib` prior to running. Note that sometimes this doesn't find all circular contigs - if it isn't working well for you, please run the original `find_circular.py`! That script is also better tested on more real world metagenomes than this one.

Usage:

```
python find_circular_fast.py contigs.fasta 150 2 contigs
```
Where 150 bp is the read length, 2 is the number of allowed mismatches, and the final argument is either 'contigs' or 'ids' to specify whether to output the circular contigs or just their ids.

Note: if 1 read length doesn't work for your assembly try 1/2 the read length - some assemblers chop off the overlap this way.


**`find_circular.py`**

Finds circular contigs in metagenome assemblies by identifying forward read overlaps at the start / end of contigs. Tested successfully on circular contigs from metagenome assemblies produced by Soapdenovo2, IDBA_UD, and SPAdes.

Requirements: (1) [Lastz](http://www.bx.psu.edu/~rsharris/lastz/) should be in your /usr/bin path, (2) BioPython.

```
usage: find_circular.py [-h] -i INPUT [-l READ_LENGTH] [-m MIN_CONTIG_SIZE]

Finds circular contigs using lastz.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input assembly filename
  -l READ_LENGTH, --read_length READ_LENGTH
                        Read length (default: 101 bp)
  -m MIN_CONTIG_SIZE, --min_contig_size MIN_CONTIG_SIZE
                        Minimum contig size to check for (default: 3 kbp)
```


**`classify.py`**

Classifies contigs based on gene size, strandedness, and intergenic region length. Works best on contigs > 20 Kb; will only attempt to classify contigs with at least 10 coding regions.

Requirements: BioPython, Scikit-learn, [Prodigal](https://github.com/hyattpd/Prodigal/wiki/installation).

```
python classify.py -i metagenome_assembly.fa
```

```
usage: classify.py [-h] -i INPUT [-m MIN_CONTIG_SIZE] [-p PRODIGAL_PATH]

Classifies viral metagenomic contigs using a random forest classifier built on
public datasets. Features are gene length, intergenic space length,
strandedness, and prodigal calling gene confidence.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input assembly filename
  -m MIN_CONTIG_SIZE, --min_contig_size MIN_CONTIG_SIZE
                        Minimum contig size to use (default: 10 kbp)
  -p PRODIGAL_PATH, --prodigal_path PRODIGAL_PATH
                        Path to prodigal (default: prodigal)
```


**`identify_host.py`**

Identifies putative bacterial and archaeal hosts for provided viral genomes by comparing tetranucleotide frequencies of viral genomes to (a) known GenBank genomes, (b) metagenomic assembled contigs with marker proteins, or (c) GenBank genomes which have been identified in the provided metagenome. Visualizes putative phage-host relationships.

Requirements: scikit-learn, hmmer3 (hmmsearch), biopython, blastp.

There are three sets of results that are calculated. The first file is `*_genbank.txt` - this approach is a naive host prediction in which the viral contigs were compared to all possible GenBank genomes, and the closest reported cellular genomes were reported. The second set of results is `*_markercontigs.txt` - this approach searches the provided entire metagenome assembly for contigs with core marker ribosomal proteins and BLASTs these proteins against a marker protein database. Viral contig tetranucleotide frequencies are then compared to the tetranucleotide frequencies of these marker host contigs, and the closest host contigs with their closest BLASTP annotation are reported for each viral contig. The third set of results is `*_markergenbank*` - this approach is similar to the 2nd, but instead of comparing viral contigs directly to metagenomic host contigs, they are compared to the full GenBank genomes of the identified bacterial / archaeal genera in the metagenomic assembly (from marker protein BLASTP results).

If no metagenome assembly is provided, then results from just the first GenBank approach will be obtained. Which method to use depends on each circumstance, but in general use the `markercontigs` results if your metagenomic assembly is sufficiently large to have many host contigs > 20 kbp.

```
python identify_host.py -i ./phage_contigs.fa -a assembled_contigs.fa -v -o ./host_identification/
```

```
usage: identify_host_markers.py [-h] -i INPUT [-a ASSEMBLY] [-v] [-b] [-hm]
                                [-o]

Predicts host(s) for a contig through tetranucleotide similarity.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Viral contig(s) fasta file
  -a ASSEMBLY, --assembly ASSEMBLY
                        Cellular metagenome assembly FASTA file (to search for
                        host marker genes).
  -v, --visualize       Visualizes NMDS of tetranucleotide frequencies for
                        host contigs and viral contigs.
  -b, --blast_path      Path to the BLASTP executable (default: blastp).
  -hm, --hmmsearch_path
                        Path to the hmmsearch executable (default: hmmsearch).
  -o, --output          Path to the directory to store output (will create if
                        does not exist).
```

**`compare_tetramers.py`**

Calculates the euclidean distance between the tetranucleotide frequencies for two given contigs or sets of contigs.

Requirements: scikit-learn.

```
python compare_tetramers.py contigs1.fa contigs2.fa
```

**`crispr_matches.py`**

Finds CRISPR loci in an assembled metagenome using minced, and searches for spacer matches to identified spacers throughout the rest of the assembled metagenome.

Requirements: BLAST+, [minced](https://github.com/ctSkennerton/minced/tree/master).

```
python crispr_matches.py -i contigs.fa
```
