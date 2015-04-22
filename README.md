tardis
======

Toolkit for Automated and Rapid DIscovery of Structural variants

Developers:

Can Alkan	calkan@gmail.com

Fereydoun Hormozdiari	hormozdiari.fereydoun@gmail.com

Emre Karakoç	ekarakoc@gmail.com

Iman Hajirasouliha	iman.ha@gmail.com

Can Koçkan	cankockan92@gmail.com

Vineet Bhakhar	vineet3692@gmail.com

Sina Jafarzadeh	jafarzadeh91@gmail.com


Building the repeats file
=========================

Download the RepeatMasker out files from the UCSC Genome Browser. For GRCh37, this file is at: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromOut.tar.gz

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

Auxiliary files
===============

GRCh37 annotations available under aux/

 * build37.dups.bed: Segmental duplication coordinates
 * build37.gaps.bed: Assembly gap coordinates
 * build37.reps.bed: RepeatMasker annotations (as described above)

Also download the reference genome from the UCSC Genome Browser. For GRCh37, this file is at: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz. Merge all FASTA files into a single file. Make sure that the same reference was used to align the reads beforehand (BAM file).


Running TARDIS
==============

	tardis -i myinput.bam --ref human_g1k_v37.fasta --gaps build37.gap.bed --reps build37.reps.bed \\
		--dups build37_dups.bed --xy --vh 

Additional parameters, helpful when debugging:

	--skip-fastq --skip-sort --skip-remap
