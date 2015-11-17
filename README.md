## Detecting Viral and Circular Contigs in Metagenomes.

**`find_circular.py`**

A python script that finds circular contigs in metagenome assemblies using two different methods.

Requirements:

1. Lastz

2. Bowtie2

**`identify_host.py`**

A python script that calculates tetranucleotide frequencies for a given contig, and returns the top ten nearest neighbors (euclidean distance) of all bacterial/archaeal genomes in GenBank. (10,000 total).

Requirements: scikit-learn.

```
python identify_post.py ./phage_contig.fa ./tetramer_database.dat
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
