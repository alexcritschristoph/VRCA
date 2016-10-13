## VIral and Circular content from metAgenomes (VICA)

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

Classifies contigs based on gene size, strandedness, and intergenic region length. Works best on contigs > 10 Kb; will only attempt to classify contigs with at least 10 coding regions.

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

Calculates tetranucleotide frequencies for a given contig, and returns the top ten nearest neighbors (euclidean distance) of all bacterial/archaeal genomes in GenBank. (10,000 total).

Requirements: scikit-learn.

```
python identify_host.py ./phage_contig.fa ./tetramer_database.dat
```

**`identify_host_markers.py`**

Calculates tetranucleotide frequencies for given viral contigs, and compares them to tetranucleotide frequencies of contigs with conserved ribosomal marker proteins from a metagenome assembly. Visualizes putative phage-host relationships.

```
python identify_host_markers.py -i ./viral_contigs.fna -a metagenome_assembly.fna -p m
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
