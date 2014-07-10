#include "tardis.h"
#define DEBUGMODE

int main( int argc, char** argv)
{
	bam_info** in_bams;
	parameters* params;
	configuration* cfg;
	int return_value;
	int i;
	int j;

	/* Load configuration file (created if it does not exist) */
	cfg = ( configuration*) malloc( sizeof( configuration));
	load_config( cfg);

	/* Set program parameters */
	init_params( &params);

	/* Parse command line arguments */	
	return_value = parse_command_line( argc, argv, params);
	if( return_value == 0)
	{
		exit( EXIT_SUCCESS);
	}
	else if( return_value != 1)
	{
		exit( return_value);
	}

	#ifdef DEBUGMODE
		print_params( params);
	#endif

	/* Read BAM files and calculate the median/avg/std of fragment sizes per library */
	in_bams = ( bam_info**) malloc( sizeof( bam_info*));
	for( i = 0; i < params->num_bams; i++)
	{
		in_bams[i] = ( bam_info*) malloc( sizeof( bam_info));  
		load_bam( in_bams[i], params->bam_file_list[i]);

		/* BAM is loaded, min/max/avg/std are calculated. Now, extract FASTQs of discordants, OEAs, and orphans */
		create_fastq( in_bams[i], params->bam_file_list[i], params);		
	}

        fprintf( stderr, "\nAll FASTQ files ready for remapping.\n");

	/* Remap with mrFAST */

	if ( params->skip_remap == 0)
	{
		return_value = remap_mrfast( params, in_bams, cfg);

		if (return_value != RETURN_SUCCESS)
			return EXIT_EXTERNAL_PROG_ERROR;
	}
	

	return EXIT_SUCCESS;
}
