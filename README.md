## VIral and Circular content from metAgenomes

**`find_circular.py`**

A python script that finds circular contigs in metagenome assemblies by identifying forward read overlaps at the start / end of contigs. Tested successfully on circular contigs from metagenome assemblies produced by Soapdenovo2, IDBA_UD, and SPAdes.

Requirements:

1. Lastz should be in your /usr/bin path.
2. BioPython should be installed.

**`classify.py`**

Classifies contigs based on gene size, strandedness, and intergenic region length. Works best on contigs > 10 Kb; will only attempt to classify contigs with at least 10 coding regions.

Requirements:

1. BioPython 
2. Scikit-learn
3. Prodigal

```
python classify.py metagenome_assembly.fa
```

**`identify_host.py`**

A python script that calculates tetranucleotide frequencies for a given contig, and returns the top ten nearest neighbors (euclidean distance) of all bacterial/archaeal genomes in GenBank. (10,000 total).

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

A python script that calculates the euclidean distance between the tetranucleotide frequencies for two given contigs or sets of contigs.

Requirements: scikit-learn.

```
python compare_tetramers.py contigs1.fa contigs2.fa
```

**`crispr_matches.py`**

A python script that finds CRISPR loci in an assembled metagenome using minced, and searches for spacer matches to identified spacers throughout the rest of the assembled metagenome.

Requirements: BLAST+, [minced](https://github.com/ctSkennerton/minced/tree/master).

```
python crispr_matches.py -i contigs.fa
```
