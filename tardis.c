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
  
  inBam->bam = gfOpen(params->bamFile, "r");
  loadBAM(inBam);
  
}
