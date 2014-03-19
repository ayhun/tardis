#include <stdio.h>
#include <stdlib.h>
#include "processbam.h"

#define SAMPLEFRAG 1000000 /*sample this many concordants to calculate avg/std*/

void loadBAM(bamInfo *inBam){
  /*
    loadBAM: Load a BAM file and parse:
    1 - Sample name (@RG SM)
    2 - Chromosome names and lengths (@SQ SN and LN)
    3 - Average fragment size; sampled from ~1M reads (SAMPLEFRAG)
    4 - Fragment size standard deviation; sampled from ~1M reads (SAMPLEFRAG)    
    
    Calculate:
    1 - min/max values. Do some sort of percentile analysis to decide std cutoff
   */

  
}


