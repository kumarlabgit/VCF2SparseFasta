# VCF2SparseFasta #

Script for extracting SNPs from a multi-sample VCF into a set of FASTA files suitable for input to MyESL program.

## Setup ##

Install bcftools: https://samtools.github.io/bcftools/howtos/install.html

Clone this repo:
	git clone https://github.com/kumarlabgit/VCF2SparseFasta

Replace line 6 in VCF2SparseFasta.py with the path to the installed bcftools executable.



## Usage ##

	cd VCF2SparseFasta
	python VCF2SparseFasta.py input_file_chr21.vcf.gz 21

By default, output will be split into FASTA files with at most 5000 positions, this can be changed with the --chunk_size parameter.

A specific range within the file can be queried by specifying --start and --end parameters.

Should accept any multi-sample VCF file that bcftools can parse, including gzip-compressed files that have been indexed.
