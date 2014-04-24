#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
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

void print_params(parameters *params)
{
	printf( "bam_file: %s\n", params->bam_file);
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

FILE* gfOpen( char* fname, char* mode)
{
	/* gfOpen: graceful file open. Try to open a file; exit if file does not exist */
	FILE* file;
	char err[500];
	file = fopen( fname, mode);
  
	if( file == NULL)
	{
		sprintf( err, "[TARDIS INPUT ERROR] Unable to open file %s in %s mode.", fname, mode[0]=='w' ? "write" : "read");
		print_error( err);
	}
	return file;
}
