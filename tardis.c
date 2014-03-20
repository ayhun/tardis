#include "tardis.h"
#define DEBUGMODE

int main(int argc, char **argv){
  parameters *params;
  int ret;
  bamInfo *inBam;

  init_params(&params);

  ret = parseCommandLine(argc, argv, params);
  if (ret == 0)
    exit(1);

  #ifdef DEBUGMODE
  print_params(params);
  #endif

  /* read BAM file and calculate avg/std of fragment sizes */

  inBam = (bamInfo *) malloc(sizeof(bamInfo));
  
  loadBAM(inBam);

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
