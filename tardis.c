#include "tardis.h"

int main( int argc, char** argv)
{
	bam_info** in_bams;
	parameters* params;
	configuration* cfg;
	int return_value;
	char username[MAX_SEQ];
	int i;
	int j;

	print_quote();
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

	if ( TARDIS_DEBUG)
        {
 	        print_params( params);
	}

	print_quote();

	/* Read BAM files and calculate the median/avg/std of fragment sizes per library */
	in_bams = ( bam_info**) malloc( sizeof( bam_info*));
	for( i = 0; i < params->num_bams; i++)
	{
		in_bams[i] = ( bam_info*) malloc( sizeof( bam_info));
		in_bams[i]->sample_name = NULL;
		load_bam( in_bams[i], params->bam_file_list[i]);

	/* BAM is loaded, min/max/avg/std are calculated. Now, extract FASTQs of discordants, OEAs, and orphans */
		if ( params->skip_fastq == 0)
	        {
		  create_fastq( in_bams[i], params->bam_file_list[i], params);		
		}
		else
		  {
		    /* TODO: check if the FASTQ files indeed exist, so it is safe to skip */
		    fprintf( stderr, "Skipping FASTQ extraction step.\n");
		}
	}
  
        fprintf( stderr, "All FASTQ files ready for remapping.\n");

	/* Remap with mrFAST */


	if ( params->skip_remap == 0)
	{
		return_value = remap_mrfast( params, in_bams, cfg);

		if (return_value != RETURN_SUCCESS)
			return EXIT_EXTERNAL_PROG_ERROR;
	}
	else
	{
	  /* TODO: check if the remapping output indeed exists, so it is safe to skip */
	        fprintf( stderr, "Skipping remapping step.\n");
	}

	print_quote();
	
	/*
	if ( params->skip_remap == 0)
	{
		return_value=divettovcf(params, in_bams);
		if (return_value != RETURN_SUCCESS){
			fprintf(stderr, "ERROR while converting DIVET files to VCF format. Exiting.\n");	
			return EXIT_EXTERNAL_PROG_ERROR;
		}
	}
	*/

	if ( params->run_vh == 1)
	{
		return_value = run_vh( params, in_bams );
		if (return_value != RETURN_SUCCESS)
			return EXIT_EXTERNAL_PROG_ERROR;
	}
	
	
	getlogin_r(username, MAX_SEQ);
	fprintf( stderr, "\n%s, before I go, I just want to tell you: you were fantastic. Absolutely fantastic. And you know what? So was I.\n", username);

	return EXIT_SUCCESS;
}
