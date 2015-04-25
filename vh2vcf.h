#ifndef __VHtoVCF
#define __VHtoVCF
// #include "tardis.h"

void do_chrom(char *thischrom, int i, int minsup, int maxsup, int mininvsup, float maxdup, float prune);
int conversion(parameters *params, bam_info ** in_bams, int i, int j);
int divettovh( parameters *params, bam_info ** in_bams);


#endif
