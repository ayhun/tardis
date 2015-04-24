#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "processfq.h"
 
void fastq_match( char* filename1, char* filename2, int num_seq, int read_length)
{
	struct read** reads1; 
	struct read** reads2;
	struct read** res;
  
	int num_batch;
	int old_batch;
	int num_matched;
	int i;
	int j;

	/*	
	FILE* f1;
	FILE* f2;
	FILE* of1;
	FILE* of2;
	*/

	gzFile  f1;
	gzFile f2;
	gzFile of1;
	gzFile of2;

	int loaded;
	char* ofilename1;
	char* ofilename2;

	ofilename1 = ( char*) malloc ( sizeof( char) * ( strlen( filename1) + 10));
	ofilename2 = ( char*) malloc ( sizeof( char) * ( strlen( filename2) + 10));
  
	sprintf( ofilename1, "%s.sorted.gz", filename1);
	sprintf( ofilename2, "%s.sorted.gz", filename2);
  
	of1 = gzopen( ofilename1, "w");
	of2 = gzopen( ofilename2, "w");

	num_batch = MEMUSE / ( read_length * 3 * 2); /* 3 comes from qname,seq,qual. 2 comes from pairs */

	if( num_batch > num_seq)
	{
		num_batch = num_seq;
	}

	alloc_reads( &reads1, num_batch);
	alloc_reads( &reads2, num_batch);

	num_matched = 0;

	f1 = gzopen( filename1, "r"); //replace with gfOpen
	f2 = gzopen( filename2, "r");

	while( num_matched < num_seq)
	{
		/* load the /1 reads */
		loaded = load_reads( f1, reads1, num_batch);
		if( loaded == 0)
		{
			/* 
			no new reads are loaded -- the array is full and it can't find more matches
			extend memory by MEMSCALE so you can load more
			*/
			old_batch = num_batch;
			num_batch = num_batch * MEMSCALE;
			if( num_batch > num_seq)
			{
				num_batch = num_seq;
			}

			realloc_reads( &reads1, old_batch, num_batch);
			realloc_reads( &reads2, old_batch, num_batch);
			loaded = load_reads( f1, reads1, num_batch);
		}    
  
		qsort( reads1, num_batch, sizeof( struct read*), fastq_qname_comp); 

		/* load the /2 reads */
		loaded = load_reads( f2, reads2, num_batch);
    
		/* search for the /2 reads within the sorted array of /1 reads */
		for( i = 0; i < num_batch; i++)
		{
			if( reads2[i]->qname != NULL && reads2[i]->empty == 0)
			{
				res = ( struct read**) bsearch( &( reads2[i]), reads1, num_batch, sizeof( struct read*), fastq_qname_comp);
				if( res != NULL)
				{
				        gzprintf( of1, "%s/1\n%s\n+\n%s\n", ( *res)->qname, ( *res)->seq, ( *res)->qual);
					gzprintf( of2, "%s/2\n%s\n+\n%s\n", reads2[i]->qname, reads2[i]->seq, reads2[i]->qual);
					( *res)->empty = 1;
					reads2[i]->empty = 1;
					num_matched++;
				}
			}
		}
	}
  
	gzclose(of1);
	gzclose(of2);
	gzclose(f1);
	gzclose(f2);

	/* 
		remove the old files (filename1, filename2)
		and rename the ofilename1 with filename1 and ofilename2 with filename2 
		so everything seems to be in-place & temp files are discarded
	*/
	
	int retval = unlink( filename1);
	if( retval == -1)
	{
		perror( "Couldn't delete tempfile1");
	}

	retval = unlink( filename2);
	if( retval == -1)
	{
		perror( "Couldn't delete tempfile2");
	}

	retval = rename( ofilename1, filename1);
	if( retval != 0) 
	{
		perror( "Failed to rename outfile1");
	}

	retval = rename( ofilename2, filename2);
	if( retval != 0) 
	{
		perror( "Failed to rename outfile2");
	}
	
	free( ofilename1);
	free( ofilename2);

	free_reads( &reads1, num_batch);
	free_reads( &reads2, num_batch);
}

static int fastq_qname_comp( const void* p1, const void* p2)
{
	struct read* a; 
	struct read* b;
	int i;
	int l1;
	int l2;

	a = *( ( struct read**) p1);
	b = *( ( struct read**) p2);

	return strcmp( a->qname, b->qname);
}

int load_reads( gzFile f1, struct read** reads, int num_batch)
{
	int i;
	int items_read;
	char qname[MAX_SEQ];
	char seq[MAX_SEQ];
	char qual[MAX_SEQ];
	char plus[MAX_SEQ];
	int cnt=0;

	for( i = 0; i < num_batch; i++)
	{
		if( reads[i]->empty == 1)
		{
			if( !gzeof( f1))
			{
			  //items_read = fscanf( f1, "@%s\n%s\n+\n%s\n", qname, seq, qual);
 			        gzgets(f1, qname, MAX_SEQ);
				if (gzeof( f1))
				  return cnt;
 			        gzgets(f1, seq, MAX_SEQ);
 			        gzgets(f1, plus, MAX_SEQ);
 			        gzgets(f1, qual, MAX_SEQ);
				/*
				items_read = sscanf( read_line, "@%s\n%s\n+\n%s\n", qname, seq, qual);
				if( items_read != 3)
				{
					perror( "Invalid number of fields read");
					}*/
				qname[strlen( qname) - 3] = 0; // get rid of /1 /2 and \n
				seq[strlen( seq) - 1] = 0; // get rid of \n
				qual[strlen( qual) - 1] = 0; // get rid of \n
				set_str( &( reads[i]->qname), qname);
				set_str( &( reads[i]->seq), seq);
				set_str( &( reads[i]->qual), qual);
				reads[i]->empty = 0;
				cnt++;
			}
		}
	}

	return cnt;
}

void alloc_reads( struct read*** reads, int num_batch)
{
	/* Allocate memory for the read data structure */
	int i;
	
	*reads = ( struct read**) malloc( sizeof( struct read*)* num_batch);

	for( i = 0; i < num_batch; i++)
	{
		( *reads)[i] = ( struct read *) malloc( sizeof( struct read));
		( *reads)[i]->qname = NULL;
		( *reads)[i]->seq = NULL;
		( *reads)[i]->qual = NULL;
		( *reads)[i]->empty = 1;
	}
}

void realloc_reads( struct read*** reads, int old_batch, int num_batch)
{
	/* Reallocate more memory for the read data structure when needed */
	int i;
	struct read** tmp_read;

	alloc_reads( &tmp_read, num_batch);

	for( i = 0; i < old_batch; i++)
	{
		set_str( &( tmp_read[i]->qname), ( *reads)[i]->qname);
		set_str( &( tmp_read[i]->seq), ( *reads)[i]->seq);
		set_str( &( tmp_read[i]->qual), ( *reads)[i]->qual);
		tmp_read[i]->empty = ( *reads)[i]->empty;
	}

	free_reads( reads, old_batch);
	*reads = tmp_read;
}

void free_reads( struct read*** reads, int num_batch)
{
	/* Free memory from the read data structure */
	int i;
	for( i = 0; i < num_batch; i++)
	{
		free( ( *reads)[i]->qname);
		free( ( *reads)[i]->seq);
		free( ( *reads)[i]->qual);
		free( ( *reads)[i]);    
	}
	free( ( *reads));
}

void create_fastq_library( struct library_properties* in_lib, char* sample_name, char* bam_path, parameters* params)
{
  /*
	FILE* fastq;
	FILE* fastq2;
	FILE* outfastq;
  */
	gzFile fastq;
	gzFile fastq2;
	gzFile outfastq;
	htsFile* bam_file;
	bam1_core_t bam_alignment_core;
	bam1_t*	bam_alignment;
	char qname[255];
	char sequence[MAX_SEQ];
	char qual[MAX_SEQ];
	char filename[255];
	char filename2[255];
	char* current_lib_name = NULL;
	int num_seq;
	int flag;
	int min;
	int max;
	int return_value;
	int i;

	/* Set FASTQ file names */
	sprintf( filename, "%s_%s_remap_1.fastq.gz", sample_name, in_lib->libname);
	sprintf( filename2, "%s_%s_remap_2.fastq.gz", sample_name, in_lib->libname);

	set_str( &( in_lib->fastq1), filename);
	set_str( &( in_lib->fastq2), filename2);

	/* if skip-fastq is set, return */
	if( params->skip_fastq)
	{
	  /* check if it is safe to skip */
		return;
	}

	/* Open FASTQ file for writing */
	//fastq = fopen( filename, "w");
	fastq = gzopen( filename, "w");
	if( !fastq)
	{
		fprintf( stderr, "Error opening the first FASTQ file\n");
		exit( 1);
	}

	/* Open the second FASTQ file for writing */	
	fastq2 = gzopen( filename2, "w");
	if( !fastq2)
	{
		fprintf( stderr, "Error opening the second FASTQ file\n");
		exit( 1);
	}

	/* Open BAM file for reading */
	bam_file = safe_hts_open( bam_path, "r");

	/* Get past the BAM header, otherwise header text gets mixed; possible htslib issue */
	bam_hdr_read( ( bam_file->fp).bgzf);

	/* Initialize the number of paired-end sequences that are discordant for the current library to zero */
	num_seq = 0;

	/* Read the BAM file alignment by alignment */
	bam_alignment = bam_init1();
	return_value = bam_read1( ( bam_file->fp).bgzf, bam_alignment);
	bam_alignment_core = bam_alignment->core;

	/* Set the read length for the current library */
	in_lib->read_length = bam_alignment_core.l_qseq;

	while( return_value != -1)
	{		
		flag = bam_alignment_core.flag;

		min = in_lib->conc_min;
		max = in_lib->conc_max;

		/* Get the library id/name for the current alignment. +1 gets rid of the leading 'Z' */
		set_str( &current_lib_name, bam_aux_get( bam_alignment, "RG") + 1);

		/* If the read is not concordant AND belongs to the current library, write it to the FASTQ file */
		if( !is_concordant( bam_alignment_core, min, max) && ( flag & BAM_FPAIRED) != 0 && strcmp( in_lib->libname, current_lib_name) == 0)
		{
			/* Increment the number of sequences written to FASTQs */			
			num_seq = num_seq + 1;
			
			/* Line 1: Read Name */
			strncpy( qname, bam_get_qname( bam_alignment), bam_alignment_core.l_qname);

			/* Read 1 goes to /1 */
			if( ( flag & BAM_FREAD1) != 0)
			{
				outfastq = fastq;
				gzprintf( outfastq, "@%s/1\n", qname);
			}
			else if( ( flag & BAM_FREAD2) != 0)
			{
				/* Read 2 goes to /2 */
				outfastq = fastq2;
				gzprintf( outfastq, "@%s/2\n", qname);
			}

			/* Line 2: Sequence */
			strncpy( sequence, bam_get_seq( bam_alignment), bam_alignment_core.l_qseq);
			sequence[bam_alignment_core.l_qseq] = '\0';

			char next_char;
			/* Read mapped to the + strand */
			if( ( flag & BAM_FREVERSE) == 0)
			{
				for( i = 0; i < strlen( sequence); i++)
				{
					next_char = base_as_char( bam_seqi( sequence, i));
					gzprintf( outfastq, "%c", next_char);
				}		       
			}
			else
			{
				/* Read mapped to the - strand */
				for( i = strlen( sequence) - 1; i >= 0; i--)
				{
					next_char = complement_char( base_as_char( bam_seqi( sequence, i)));
					gzprintf( outfastq, "%c", next_char);
				}		       
			}
			
			/* Line 3: "+" */
			gzprintf( outfastq, "\n+\n");

			/* Line 4: Quality String */
			strncpy( qual, bam_get_qual( bam_alignment), bam_alignment_core.l_qseq);
			qual[bam_alignment_core.l_qseq] = '\0';

			/* Convert the quality value to ASCII */
			qual_to_ascii( qual);
		
			/* If the read is mapped to the reverse strand, reverse the quality string */
			if( bam_is_rev( bam_alignment))
			{
			  	reverse_string( qual);
			}
			gzprintf( outfastq, "%s\n", qual);
		}

		/* Read next alignment */
		return_value = bam_read1( ( bam_file->fp).bgzf, bam_alignment);
		bam_alignment_core = bam_alignment->core;
	}

	/* We should divide num_seq by 2 as the reads are paired and half of the total will be written to each FASTQ file */
	num_seq = num_seq / 2;

	/* Set the number of sequences that are written to each FASTQ file for the current library */
	in_lib->num_sequences = num_seq;
		
	/* Close the BAM file */
	return_value = hts_close( bam_file);
	if( return_value != 0)
	{
		fprintf( stderr, "Error closing BAM file\n");
		exit( 1);
	}

	/* Close the FASTQ file */
	return_value = gzclose( fastq);
	if( return_value != 0)
	{
		fprintf( stderr, "Error closing the first FASTQ file\n");
	}

	/* Close the second FASTQ file */
	return_value = gzclose( fastq2);
	if( return_value != 0)
	{
		fprintf( stderr, "Error closing the second FASTQ file\n");
	}

	/* Free memory */
	free( current_lib_name);

	if( !( params->skip_sort))
        {
	        fprintf( stderr, "Sorting FASTQ files for library: %s.\n", in_lib->libname);
	        fastq_match( in_lib->fastq1, in_lib->fastq2, in_lib->num_sequences, in_lib->read_length);
	}	    
	
}
