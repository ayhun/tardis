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
	char cmdline[2048];

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
		for( i = 0; i < params->num_bams; i++)
		{
			for ( j = 0; j < in_bams[i]->num_libraries; j++)
		        {
				fprintf( stderr, "\nRemapping:\n\t\tSample: %s\n\t\tLibrary: %s\n", in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname);
				sprintf( cmdline, "%s --search %s --pe --min %d --max %d --sample %s --lib %s --rg %s --seq1 %s --seq2 %s -o %s-%s.sam -u %s-%s.unmapped.fastq", 
					cfg->path_mrfast, params->ref_genome, in_bams[i]->libraries[j]->conc_min, in_bams[i]->libraries[j]->conc_max,
					in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname, in_bams[i]->libraries[j]->libname, 
					in_bams[i]->libraries[j]->fastq1, in_bams[i]->libraries[j]->fastq2,
					in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname,
					in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname);
				return_value = system( cmdline);
				
				if (WIFSIGNALED(return_value) && (WTERMSIG(return_value) == SIGINT || WTERMSIG(return_value) == SIGQUIT))
				{
					fprintf( stderr, "mrFAST remapping failed.\n");
					return EXIT_EXTERNAL_PROG_ERROR;
				}

				/* to be implemented.
					params->threads will be used for multithreading option of mrFAST
				*/		       
			}
		}  
	}


	return EXIT_SUCCESS;
}
