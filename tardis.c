#include "tardis.h"
#define DEBUGMODE

int main( int argc, char** argv)
{
	bam_info* in_bam;
	parameters* params;
	int return_value;

	/* Set program parameters */
	init_params( &params);

	/* Parse command line arguments */	
	return_value = parse_command_line( argc, argv, params);
	if( return_value == 0)
	{
		exit( 1);
	}

	#ifdef DEBUGMODE
		print_params( params);
	#endif

	/* Read BAM file and calculate the median/avg/std of fragment sizes */
	in_bam = ( bam_info*) malloc( sizeof( bam_info));
  
	load_bam( in_bam, params->bam_file);

	/*#ifdef DEBUGMODE
		print_bam( in_bam);
	#endif*/
	

  /*
    BAM is loaded, min/max/avg/std are calculated.
    Now, extract FASTQs of discordants, OEAs, and orphans
  */
  
  /*  to be implemented.  We need to decide on the FASTQ name convention.
    createFastqs(inBam);
  */

  /* 
     remap with mrFAST
  */
  
  /*  to be implemented.
      pass config (mrfast path) and FASTQ names. 
      params->threads will be used for multithreading option of mrFAST
      remap(params, ...); 
   */
}
