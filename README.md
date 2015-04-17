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

Concetanate all out files and get rid of unneeded information:

	cat `find rmasker -name *.out` \
		| grep -v "matching\|begin"\
		| awk '{OFS="\t"; if ($1 ~ /./) print $5,$6-1,$7,$9,$10,$11}'\
		| sed s/chr// | sort -k 1,1 -k 2,2n > build37.reps.bed

Remove unnecessary files:

	rm chromOut.tar.gz
	rm -fr rmasker

