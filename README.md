tardis
======

Toolkit for Automated and Rapid DIscovery of Structural variants

Requirements
============

 * zlib   (http://www.zlib.net)
 * mrfast (https://github.com/BilkentCompGen/mrfast)
 * htslib (included as submodule; http://htslib.org/)

Fetching tardis
===============

	git clone https://github.com/calkan/tardis.git --recursive

Compilation
===========

Type:

	make libs
	make
	cp tardis /path/to/your/favorite/binaries


Auxiliary files
===============

GRCh37 annotations available under aux/

 * build37.dups.bed: Segmental duplication coordinates.
 * build37.gaps.bed: Assembly gap coordinates.
 * build37.reps.bed: RepeatMasker annotations (as described below). This file is provided as compressed. Unzip it before use.

Also download the reference genome from the UCSC Genome Browser. For GRCh37, this file is at: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz. Merge all FASTA files into a single file. Make sure that the same reference was used to align the reads beforehand (BAM file).

** Reference genome and its annotations should use the EXACT same names for the chromosomes. The example provided below use the 1000 Genomes Project naming convention. **

Building the repeats file
=========================

Download the RepeatMasker out files from the UCSC Genome Browser. For GRCh37 (hg19), this file is at: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromOut.tar.gz

Create a temp directory to make things easier to cleanup:

	mkdir rmasker
	tar zxf chromOut.tar.gz -C rmasker

Concetanate all out files and get rid of unneeded information. Make sure that the chromosome names in the reference genome and the repeats file matches.

	cat `find rmasker -name *.out` \
		| grep -v "matching\|begin"\
		| awk '{OFS="\t"; if ($1 ~ /./) print $5,$6-1,$7,$9,$10,$11}'\
		| sed s/chr// | sed s/Un_// | sed s/_random// | sed s/gl/GL/ | grep -v hap \
		| sed s/\.\._// | sed s/._// \
		| awk '{OFS="\t"; if ($1=="M") $1="MT"; if ($1 ~ /GL/) $1=$1".1"; print $0}'\
		| sort -k 1,1 -k 2,2n > build37.reps.bed

Remove unnecessary files:

	rm chromOut.tar.gz
	rm -fr rmasker


Running tardis
==============

	tardis -i myinput.bam --ref human_g1k_v37.fasta --gaps build37.gap.bed \
		--reps build37.reps.bed --dups build37_dups.bed --xy --vh \
		--out myoutput

Additional parameters, helpful when debugging:

	--skip-fastq --skip-sort --skip-remap

All parameters
==============

	--bamlist   [bamlist file] : A text file that lists input BAM files one file per line.
	--input [BAM files]        : Input files in sorted and indexed BAM format. You can pass multiple BAMs using multiple --input parameters.
	--out   [output prefix]    : Prefix for the output file names.
	--ref   [reference genome] : Reference genome in FASTA format.
	--gaps  [gaps file]        : Assembly gap coordinates in BED3 format.
	--dups  [dups file]        : Segmental duplication coordinates in BED3 format.
	--reps  [reps file]        : RepeatMasker annotation coordinates in BED6 format. See manual for details.
	--mei   ["Alu:L1Hs:SVA"]   : List of mobile element names.
	--xx                       : Sample is male.
	--xy                       : Sample is female.
	--vh                       : Run VariationHunter/CommonLAW (read pair + read depth).
	--ns                       : Run NovelSeq (read pair + assembly).
	--sr                       : Run SPLITREAD (split read).
	--all                      : Run all three algorithms above [DEFAULT].
	--skip-fastq               : Skip FASTQ dump for discordants. Use this only if you are regenerating the calls.
	--skip-sort                : Skip FASTQ sort for discordants. Use this only if you are regenerating the calls.
	--skip-remap               : Skip FASTQ remapping for discordants. Use this only if you are regenerating the calls.
	--version                  : Print version and exit.
	--help                     : Print this help screen and exit.

