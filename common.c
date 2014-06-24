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
	( *params)->bam_file = NULL;
	( *params)->gaps = NULL;
	( *params)->mei = NULL;
	( *params)->sample_gender = MALE;
	( *params)->run_vh = 0; 
	( *params)->run_ns = 0;
	( *params)->run_sr = 0;
	( *params)->threads = 1;
}

void print_params( parameters *params)
{
	printf( "bam_file: %s\n", params->bam_file);
	printf( "ref_genome: %s\n", params->ref_genome);
	printf( "reps: %s\n", params->reps);
	printf( "dups: %s\n", params->dups);
	printf( "gaps: %s\n", params->gaps);
	printf( "mei: %s\n", params->mei);
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
	char library_buffer[1024];
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

			set_str( &( in_bam->( libraries[i])->libname), library_name_buffer);
			i = i + 1;
		}

		p = strtok( NULL, "\n");
	}
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
