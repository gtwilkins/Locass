# Locass

## Overview
Locass is an application for the assembly of targeted genomic loci and their genomic neighbourhood. From a dataset unassembled sequence reads, locass will find and assemble genomic regions that are similar to a user's query sequence. Locass will extend the assembled loci into flanking genomic regions as desired. Locass can be run with minimal memory requirements even with large eukaryotic genomes.

## Limitations
Locass is designed for use with short read sequence data, namely those produced by Illumina. Data must include at least one paired read library.

## Requirements
* gcc

## Installation
The install directory can be specified with the following command (if this omitted, Locass is installed to /usr/local/bin/):

	./configure --prefix=/path/to/directory/

Locass can be installed with the following commands:

	make
	sudo make install

## Use
Locass assembly requires three steps. The first two steps are preprocessing steps that must be performed once for a given read dataset. Once preprocessing is complete, assembly can be repeated without the first two steps.

### Step 1: Read indexing
The sequence read dataset must be indexed to allow for rapid assembly with the `index` command.

The input for `index` should be a single file containing a list of input read libraries, where each line details a library, indicating whether the reads are paired and their filename(s). Paired libraries can be either separated into two separate files or inteleaved in a single file. A given library is thus detailed in one of three forms as follows:

	paired /path/to/separated_read_file-1.fastq /path/to/separated_read_file-2.fastq
	paired /path/to/interleaved_read_file.fa
	single /path/to/unpaired_read_file

An example `index` command would be:

	locass index -i /path/to/sequence_file_list -p /path/to/project_prefix

### Step 2: Library calibration
The orientation and insert size of paired read libraries, as well as overall sequencing coverage must be detected before assembly with the `calibrate` command.

An example `calibration` command would be:

	locass calibration -p /path/to/project_prefix

### Step 3: Assembly
Targeted locus assembly can be performed once the sequence read dataset has been indexed and calibrated with the `assemble` command.

The input for `assemble should be a single file containing one or more query sequences, either the form of a fasta file or a simple list of sequences, one per line. Query sequences do not need to be long (approximately the length of a sequence read is sufficient), as assembly will extend beyond the query.

An example `assemble` command would be:

	locass assemble -i /path/to/seed_file.fa -p /path/to/project_prefix -o /path/to/output_file.fa

Additional assembly options can be viewed with:

	locass assemble -h
