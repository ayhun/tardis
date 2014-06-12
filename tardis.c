#include "tardis.h"
#define DEBUGMODE


int main( int argc, char** argv)
{
	bam_info* in_bam;
	parameters *params;
	MainOptions mainOptions;
	configuration* cfg;
	int return_value;

	/* Set program parameters */
	init_params( &params);

	/* Parse command line arguments */	
	return_value = parse_command_line( argc, argv, params, &mainOptions);
	if( return_value == 0)
	{
		exit( 1);
	}
	printf("%d\n", return_value);

	// #ifdef DEBUGMODE
		print_params( params);
	// #endif
	printf(" \n"); //temp: @VINEET
	/* Load configuration file */
	cfg = ( configuration*) malloc( sizeof( configuration));
	load_config( cfg);

	/* Read BAM file and calculate the median/avg/std of fragment sizes */
	in_bam = ( bam_info*) malloc( sizeof( bam_info));  
	load_bam( in_bam, params->bam_file);
	printf(" \n"); //temp: @VINEET	
	/* BAM is loaded, min/max/avg/std are calculated. Now, extract FASTQs of discordants, OEAs, and orphans */
	create_fastq( in_bam, params);
	printf(" \n"); //temp: @VINEET
	/* sort FASTQ files to match /1 and /2 reads; unless skip-sort is invoked
	if( !( params->skip_sort))
	{
		sortFastqs( in_bam, params);
	}
	*/
	/* Remap with mrFAST */
	run(mainOptions.libFileAdrs, mainOptions.chroFileName, mainOptions.gapFileName, mainOptions.repeatFileName, mainOptions.initializeFileName, mainOptions.prunProb, mainOptions.svSup, mainOptions.outputFile, mainOptions.outputRead, mainOptions.overMapLimit);
		
  	// vhprocessing(argc, argv,return_value);
  	printf(" \n"); //temp: @VINEET
	/* to be implemented.
	pass config (mrfast path) and FASTQ names. 
	params->threads will be used for multithreading option of mrFAST
	remap(params, ...); 
	*/
}
