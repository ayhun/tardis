#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>

/* tardis headers */
#include "processbam.h"

/* Sample this many fragments to calculate avg/median/std per library */
#define SAMPLEFRAG 1000000 

/* Maximum sequence/quality length */
#define MAX_SEQ 1000

void load_bam( bam_info* in_bam, char* path)
{
	/* Variables */
	htsFile* bam_file;
	bam_hdr_t* bam_header;
	bam1_core_t bam_alignment_core;
	bam1_t*	bam_alignment;
	int** fragment_size;
	int** second_pass_fragments;
	int* fragments_sampled;
	int* second_test_pass;
	int* fragment_size_total;
	float* variance;
	char* library_name;
	int library_index;
	int diff;
	int return_value;
	int i;
	int j;

	/* Open the BAM file for reading. htslib automatically detects the format
		of the file, so appending "b" after "r" in mode is redundant. */
	bam_file = safe_hts_open( path, "r");

	/* Read in BAM header information */
	bam_header = bam_hdr_read( ( bam_file->fp).bgzf);

	/* Store the number of reference sequences */
	in_bam->num_chrom = bam_header->n_targets;

	/* Allocate memory for reference sequence lengths */
	in_bam->chrom_lengths = ( int*) malloc( in_bam->num_chrom * sizeof( int));

	in_bam->chrom_names = ( char**) malloc( in_bam->num_chrom * sizeof( char*));
	
	/* Store chromosome lengths, allocate memory for reference sequence names,
	 and store the names as well */
	for( i = 0; i < in_bam->num_chrom; i++)
	{
		set_str( ( &( in_bam->chrom_names)[i]), ( bam_header->target_name)[i]);
	}
	
	/* Extract the Sample Name from the header text */
	get_sample_name( in_bam, bam_header->text);

	/* Extract the number of libraries within the BAM file */
	get_library_count( in_bam, bam_header->text);

	/* Initialize the structures for library properties */
	in_bam->libraries = ( struct library_properties**) malloc( in_bam->num_libraries * sizeof( struct library_properties*));
	for( i = 0; i < in_bam->num_libraries; i++)
	{
		( in_bam->libraries)[i] = ( struct library_properties*) malloc( sizeof( struct library_properties));
	}

	/* Extract the ids/names for the libraries. A single Sample with multiple 
	 possible libraries are assumed for each BAM file */
	get_library_names( in_bam, bam_header->text);
 
	/* For SAMPLEFRAG number of alignments, store the template length field.
	 Performed for each different library */
	fragment_size = ( int**) malloc( in_bam->num_libraries * sizeof( int*));
	for( i = 0; i < in_bam->num_libraries; i++)
	{
		fragment_size[i] = ( int*) malloc( SAMPLEFRAG * sizeof( int));
	}

	/* Initial read */	
	bam_alignment = bam_init1();
	return_value = bam_read1( ( bam_file->fp).bgzf, bam_alignment);

	/* The remaining reads */
	fragments_sampled = ( int*) malloc( in_bam->num_libraries * sizeof( int));
	for( i = 0; i < in_bam->num_libraries; i++)
	{
		fragments_sampled[i] = 0;
	}

	while( return_value != -1 && !sufficient_fragments_sampled( fragments_sampled, in_bam->num_libraries))
	{
		bam_alignment_core = bam_alignment->core;

		if( bam_alignment_core.isize > 0 && !bam_is_rev( bam_alignment) && bam_is_mrev( bam_alignment))
		{
			set_str( &library_name, bam_aux_get( bam_alignment, "RG"));
			library_index = find_library_index( in_bam, library_name);

			/* Sample SAMPLEFRAG number of alignments for each library at most */
			if( fragments_sampled[library_index] < SAMPLEFRAG)
			{
				fragment_size[library_index][fragments_sampled[library_index]] = bam_alignment_core.isize;
				fragments_sampled[library_index] = fragments_sampled[library_index] + 1;
			}
		}

		/* Read next alignment */
		return_value = bam_read1( ( bam_file->fp).bgzf, bam_alignment);
	}

	/* Now we have SAMPLEFRAG number of fragment sizes which are positive and pass the flag conditions */
	for( i = 0; i < in_bam->num_libraries; i++)
	{
		/* Sort the fragment sizes */
		qsort( fragment_size[i], SAMPLEFRAG, sizeof( int), compare_size_int);
		
		/* Get the medians */
		( in_bam->libraries)[i]->frag_med = fragment_size[i][( SAMPLEFRAG / 2) - 1];
	}
	
	/* Find the fragment sizes which pass the second test, and will contribute to the avg and std */
	second_pass_fragments = ( int**) malloc( in_bam->num_libraries * sizeof( int*));
	for( i = 0; i < in_bam->num_libraries; i++)
	{
		second_pass_fragments[i] = ( int*) malloc( SAMPLEFRAG * sizeof( int));
	}

	second_test_pass = ( int*) malloc( in_bam->num_libraries * sizeof( int));
	fragment_size_total = ( int*) malloc( in_bam->num_libraries * sizeof( int));

	for( i = 0; i < in_bam->num_libraries; i++)
	{
		second_test_pass[i] = 0;
		fragment_size_total[i] = 0;
	}

	for( i = 0; i < in_bam->num_libraries; i++)
	{
		for( j = 0; j < SAMPLEFRAG; j++)
		{
			if( fragment_size[i][j] <= 2 * ( in_bam->libraries)[i]->frag_med)
			{
				fragment_size_total[i] = fragment_size_total[i] + fragment_size[i][j];
				second_pass_fragments[i][j] = fragment_size[i][j];
				second_test_pass[i] = second_test_pass[i] + 1;
			}
		}
	}

	for( i = 0; i < in_bam->num_libraries; i++)
	{
		/* Compute the averages */
		( in_bam->libraries)[i]->frag_avg = ( float) fragment_size_total[i] / ( float) second_test_pass[i];
	}

	/* Compute the variance and std */
	variance = ( float*) malloc( in_bam->num_libraries * sizeof( float));
	for( i = 0; i < in_bam->num_libraries; i++)
	{
		variance[i] = 0;
	}

	for( i = 0; i < in_bam->num_libraries; i++)
	{
		for( j = 0; j < second_test_pass[i]; j++)
		{
			diff = ( second_pass_fragments[i][j] - ( in_bam->libraries)[i]->frag_avg);
			variance[i] = variance[i] +  diff * diff;
		}
	}

	for( i = 0; i < in_bam->num_libraries; i++)
	{
		variance[i] = ( float) variance[i] / ( float) second_test_pass[i];
		( in_bam->libraries)[i]->frag_std = sqrt( variance[i]);
	}
	
	/* Close the BAM file */
	return_value = hts_close( bam_file);
	if( return_value != 0)
	{
		fprintf( stderr, "Error closing BAM file\n");
		exit( 1);
	}

	/* Free Memory */
	for( i = 0; i < in_bam->num_libraries; i++)
	{
		free( fragment_size[i]);
		free( second_pass_fragments[i]);
	}
	free( fragment_size);
	free( second_pass_fragments);	
	free( fragments_sampled);
	free( second_test_pass);
	free( fragment_size_total);
	free( variance);
	free( library_name);
}

void get_sample_name( bam_info* in_bam, char* header_text)
{
	/* Delimit the BAM header text with tabs and newlines */
	char* p = strtok( header_text, "\t\n");
	char sample_name_buffer[1024];

	while( p != NULL)
	{
		/* If the current token has "SM" as the first two characters,
			we have found our Sample Name */
		if( p[0] == 'S' && p[1] == 'M')
		{
			/* Get the Sample Name */
			strncpy( sample_name_buffer, p + 3, strlen( p) - 3);

			/* Add the NULL terminator */
			sample_name_buffer[strlen( p) - 3] = '\0';

			/* Exit loop */
			break;
		}
		p = strtok( NULL, "\t\n");
	}

	set_str( &( in_bam->sample_name), sample_name_buffer);
}

void get_library_count( bam_info* in_bam, char* header_text)
{
	int number_of_libraries = 0;

	/* Delimit the BAM header text with newlines */
	char* p = strtok( header_text, "\n");
	while( p != NULL)
	{
		/* If the current token (which is a line) has "RG" as the second and third characters,
		 we have found a new library */
		if( p[1] == 'R' && p[2] == 'G')
		{
			number_of_libraries = number_of_libraries + 1;
		}
	}

	in_bam->num_libraries = number_of_libraries;
}

void get_library_names( bam_info* in_bam, char* header_text)
{
	char line_buffer[1024];
	char library_name_buffer[1024];
	int i;

	/* Delimit the BAM header text with newlines */
	i = 0;
	char* p = strtok( header_text, "\n");
	while( p != NULL)
	{
		/* If the current token (which is a line) has "RG" as the second and third characters,
			we copy the current line to a new char* and tokenize it to get the id/name */
		if( p[1] == 'R' && p[2] == 'G')
		{
			/* Copy the current line to the line buffer */
			strncpy( line_buffer, p, strlen( p));
			line_buffer[strlen( p)] = '\0';
			
			/* Tokenize the line buffer by tabs */
			char* q = strtok( line_buffer, "\t");
			while( q != NULL)
			{
				if( q[0] == 'I' && p[1] == 'D')
				{
					/* Get the Library Name */
					strncpy( library_name_buffer, q + 3, strlen( q) - 3);

					/* Add the NULL terminator */
					library_name_buffer[strlen( q) - 3] = '\0';

					/* Exit loop */
					break;
				}

				q = strtok( NULL, "\t");
			}

			set_str( &( ( in_bam->libraries)[i]->libname), library_name_buffer);
			i = i + 1;
		}

		p = strtok( NULL, "\n");
	}
}

int find_library_index( bam_info* in_bam, char* library_name)
{
	int library_index;
	int i;

	for( i = 0; i < in_bam->num_libraries; i++)
	{
		if( strcmp( ( in_bam->libraries)[i]->libname, library_name) == 0)
		{
			return i;
		}
	}

	return -1;
}
	
int sufficient_fragments_sampled( int* fragments_sampled, int num_libraries)
{
	int i;
	for( i = 0; i < num_libraries; i++)
	{
		if( fragments_sampled[i] != SAMPLEFRAG)
		{
			return 0;
		}
	}
	
	return 1;
}
		
void print_bam( bam_info* in_bam)
{
	printf( "Number of Chromosomes: %d\n", in_bam->num_chrom);

	int i;
	for( i = 0; i < in_bam->num_chrom; i++)
	{
		printf( "Chromosome Name: %s\n", ( in_bam->chrom_names)[i]);
		printf( "Length of the Chromosome: %d\n", ( in_bam->chrom_lengths)[i]);
	}
}

void create_fastq( bam_info* in_bam, parameters *params)
{
	FILE* fastq;
	FILE* fastq2;
	FILE *outfastq;
	htsFile* bam_file;
	bam1_core_t bam_alignment_core;
	bam1_t*	bam_alignment;
	char qname[255];
	char sequence[MAX_SEQ];
	char qual[MAX_SEQ];
	char filename[255];
	char filename2[255];
	int flag;
	int min;
	int max;
	int return_value;
	int i;
	char *path = params->bam_file;

	/* Set FASTQ file names */
	sprintf( filename, "%s_remap_1.fastq", in_bam->sample_name);
	sprintf( filename2, "%s_remap_2.fastq", in_bam->sample_name);

	set_str( &( in_bam->fastq1), filename);
	set_str( &( in_bam->fastq2), filename2);

	/* if skip-fastq is set, return */
	if( params->skip_fastq)
	{
		return;
	}

	/* Open FASTQ file for writing */
	fastq = fopen( filename, "w");
	if( !fastq)
	{
		fprintf( stderr, "Error opening the first FASTQ file\n");
		exit( 1);
	}

	/* Open the second FASTQ file for writing */	
	fastq2 = fopen( filename2, "w");
	if( !fastq2)
	{
		fprintf( stderr, "Error opening the second FASTQ file\n");
		exit( 1);
	}

	/* Open BAM file for reading */
	bam_file = hts_open( path, "r");
	if( !bam_file)
	{
		fprintf( stderr, "Error opening BAM file\n");
		exit( 1);
	}

	/* Get past the BAM header, otherwise header text gets mixed; possible htslib issue */
	bam_hdr_read( ( bam_file->fp).bgzf);

	/* For all reads within the BAM file, check if concordant */
	bam_alignment = bam_init1();
	return_value = bam_read1( ( bam_file->fp).bgzf, bam_alignment);

	while( return_value != -1)
	{
		bam_alignment_core = bam_alignment->core;
		flag = bam_alignment_core.flag;


		min = in_bam->frag_avg - ( 4 * in_bam->frag_std);
		max = in_bam->frag_avg + ( 4 * in_bam->frag_std);
		
		if( min < 0)
		{
			min = 0;
		}

		/* If the read is not concordant, write it to the FASTQ file */
		if( !is_concordant( bam_alignment_core, min, max) && ( flag & BAM_FPAIRED) != 0)
		{
			/* Line 1: Read Name */
			strncpy( qname, bam_get_qname( bam_alignment), bam_alignment_core.l_qname);

			/* Read 1 goes to /1 */
			if( ( flag & BAM_FREAD1) != 0)
			{
				outfastq = fastq;
				fprintf( outfastq, "@%s/1\n", qname);
			}
			else if( ( flag & BAM_FREAD2) != 0)
			{
				/* Read 2 goes to /2 */
				outfastq = fastq2;
				fprintf( outfastq, "@%s/2\n", qname);
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
					fprintf( outfastq, "%c", next_char);
				}		       
			}
			else
			{
				/* Read mapped to the - strand */
				for( i = strlen( sequence) - 1; i >= 0; i--)
				{
					next_char = complement_char( base_as_char( bam_seqi( sequence, i)));
					fprintf( outfastq, "%c", next_char);
				}		       
			}
			
			/* Line 3: "+" */
			fprintf( outfastq, "\n+\n");

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
			fprintf( outfastq, "%s\n", qual);
		}

		/* Read next alignment */
		return_value = bam_read1( ( bam_file->fp).bgzf, bam_alignment);
	}
		
	/* Close the BAM file */
	return_value = hts_close( bam_file);
	if( return_value != 0)
	{
		fprintf( stderr, "Error closing BAM file\n");
		exit( 1);
	}

	/* Close the FASTQ file */
	return_value = fclose( fastq);
	if( return_value != 0)
	{
		fprintf( stderr, "Error closing the first FASTQ file\n");
	}

	/* Close the second FASTQ file */
	return_value = fclose( fastq2);
	if( return_value != 0)
	{
		fprintf( stderr, "Error closing the second FASTQ file\n");
	}
}
