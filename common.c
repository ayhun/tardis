#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

/* htslib headers */
#include <htslib/sam.h>
#include <htslib/hts.h>

#include "common.h"

void init_params( parameters** params)
{
	/* initialize parameters */
	*params = ( parameters*) malloc( sizeof( parameters));
	( *params)->ref_genome = NULL;
	( *params)->reps = NULL;
	( *params)->dups = NULL;
	( *params)->bam_files = NULL;
	( *params)->bam_list_path = NULL;
	( *params)->bam_file_list = NULL;
	( *params)->gaps = NULL;
	( *params)->mei = NULL;
	( *params)->sample_gender = MALE;
	( *params)->run_vh = 0; 
	( *params)->run_ns = 0;
	( *params)->run_sr = 0;
	( *params)->threads = 1;
	( *params)->num_bams = 0;
}

void print_params( parameters* params)
{
	int i;

	printf( "Number of bam files: %d\n", params->num_bams);
	for( i = 0; i < params->num_bams; i++)
	{
		printf( "BAM input: %s\n", params->bam_file_list[i]);
	}
	printf( "ref_genome: %s\n", params->ref_genome);
	printf( "reps: %s\n", params->reps);
	printf( "dups: %s\n", params->dups);
	printf( "gaps: %s\n", params->gaps);
	printf( "mei: %s\n", params->mei);
}

void print_error( char* msg)
{
	/* print error message and exit */
	fprintf( stderr, "\n%s\n", msg);
	fprintf( stderr, "Invoke parameter -h for help.\n");
	exit( 1);
}

FILE* safe_fopen( char* path, char* mode)
{
	/* Safe file open. Try to open a file; exit if file does not exist */
	FILE* file;
	char err[500];

	file = fopen( path, mode);  
	if( !file)
	{
		sprintf( err, "[TARDIS INPUT ERROR] Unable to open file %s in %s mode.", path, mode[0]=='w' ? "write" : "read");
		print_error( err);
	}
	return file;
}

htsFile* safe_hts_open( char* path, char* mode)
{
	htsFile* bam_file;
	char err[500];
	
	bam_file = hts_open( path, mode);
	if( !bam_file)
	{
		sprintf( err, "[TARDIS INPUT ERROR] Unable to open file %s in %s mode.", path, mode[0]=='w' ? "write" : "read");
		print_error( err);
	}

	return bam_file;
}

int is_concordant( bam1_core_t bam_alignment_core, int min, int max)
{
	int flag = bam_alignment_core.flag;

	if( ( flag & BAM_FPROPER_PAIR) == 0) 
	{
		/* Not proper pair */
		return 0;
	}

	if( ( flag & BAM_FUNMAP) != 0)  // c.a.
	{
		/* Read unmapped; Orphan or OEA */
		return 0;
	}

	if( ( flag & BAM_FMUNMAP) != 0) // c.a.
	{
		/* Mate unmapped; Orphan or OEA */
		return 0;
	}

	if( ( flag & BAM_FREVERSE) != 0 && ( flag & BAM_FMREVERSE) != 0)
	{
		/* -- orientation = inversion */
		return 0;
	}

	if( ( flag & BAM_FREVERSE) == 0 && ( flag & BAM_FMREVERSE) == 0)
	{
		/* ++ orientation = inversion */
		return 0;
	}

	if( bam_alignment_core.tid != bam_alignment_core.mtid) 
	{
		/* On different chromosomes */
		return 0;
	}

	if( bam_alignment_core.pos <= bam_alignment_core.mpos) // c.a.
	{
		/* Read is placed BEFORE its mate */
		if( ( flag & BAM_FREVERSE) != 0 && ( flag & BAM_FMREVERSE) == 0)
		{
			/* -+ orientation = tandem duplication */
			return 0;
		}
	}
	else
	{
		/* Read is placed AFTER its mate */
		if( ( flag & BAM_FREVERSE) == 0 && ( flag & BAM_FMREVERSE) != 0)
		{
			/* +- orientation = tandem duplication */
			return 0;
		}
	}

	/* Passed all of the above. proper pair, both mapped, in +- orientation. Now check the isize */
	if( abs(bam_alignment_core.isize) < min || abs(bam_alignment_core.isize) > max) // c.a.
	{
		/* Deletion or Insertion */
		return 0;
	}

	/* All passed. Read is concordant */
	return 1;
}

/* Decode 4-bit encoded bases to their corresponding characters */
char base_as_char( int base_as_int)
{
	if( base_as_int == 1)
	{
		return 'A';
	}
	else if( base_as_int == 2)
	{
		return 'C';
	}
	else if( base_as_int == 4)
	{
		return 'G';
	}
	else if( base_as_int == 8)
	{
		return 'T';
	}
	else if( base_as_int == 15)
	{
		return 'N';
	}
}

/* Return the complement of a base */
char complement_char( char base)
{
	switch( base)
	{
		case 'A': 
			return 'T';
			break;
		case 'C': 
			return 'G';
			break;
		case 'G': 
			return 'C';
			break;
		case 'T': 
			return 'A';
			break;
		default: 
			return 'N';
			break;
	}
	return 'X';
}

/* Add 33 to the interger value of the qual characters to convert them to ASCII */
void qual_to_ascii( char* qual)
{
	int i;
	for( i = 0; i < strlen( qual); i++)
	{
		qual[i] = qual[i] + 33;
	}
}

/* Even safer than strncpy as it dynamically allocates space for the string if
 there hasn't been already */
void set_str( char** target, char* source)
{
	if( *target != NULL)
	{
		free( ( *target));
	}

	( *target) = ( char*) malloc( sizeof( char) * ( strlen( source) + 1));
	strncpy( ( *target), source, ( strlen( source) + 1));
}

/* Reverse a given string */
void reverse_string( char* str)
{
	int i;
	char swap;
	int len = strlen( str);

	for( i = 0; i < len / 2; i++)
	{
		swap = str[i];
		str[i] = str[len - i - 1];
		str[len - i - 1] = swap;
	}
}

int compare_size_int( const void* p, const void* q)
{
    int i = *( const int*) p;
    int j = *( const int*) q;

	if( i < j)
	{
		return -1;
	}
	else if( i == j)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}
