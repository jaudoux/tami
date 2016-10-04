# TaMI : Targeted Mutation Identification

TaMI, is a k-mer base variant caller that directly scan raw FASTQ files and looks for
variants without alignement by using a dictionnary of mutated k-mers constructed for specific targets.

TaMI looks for all possible single nucleotide mutations (1-nt indels as well) only on the targeted regions.
TaMI is very fast, a whole human exome can be scan within a minute in a single thread.
TaMI memory load depend on the size of the indexed target regions.

# Installation

- Clone tami repositiory
- run `make` to compile tami
- place the `tami` binary somewhere accessible from your `$PATH`.

# Usage

First you need to create a TAM files, wich is a dictionnary of mutated k-mers composed of
all possible mutation in the targets (bed file).

`tami build -r reference.fasta -o file.tam file.bed`

Then you need to scan FASTQ files to look for mutations indexed in the TAM.

`tami scan file.tam reads_1.fastq.gz reads_2.fastq.gz > output.vcf`

# Authors

Jérôme Audoux, Alexandre Soriano.
