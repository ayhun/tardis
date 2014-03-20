#ifndef __PROCESSBAM
#define __PROCESSBAM

typedef struct _bamInfo{
  FILE *bam; /* file pointer to the BAM file; as set in parameters->bamFile */
  /* FILE * above will be replaced by the pointer type for BAMs as defined in HTSLib */
  float fragAvg; /* average fragment size */
  float fragStd; /* fragment size standard deviation */
  int concMin; /* min cutoff for concordants */
  int concMax; /* max cutoff for concordants */
  char *sampleName; /* name of the sample, parsed from SM in the BAM header */
  /* need something for chromosome names and lengths, parse from BAM header */
  int numChrom; /* number of chromosomes */
  char **chromNames; /* names of the chromosomes */
  int *chromLengths; /* lengths of the chromosomes */
} bamInfo;


void loadBAM(bamInfo *);

#endif
