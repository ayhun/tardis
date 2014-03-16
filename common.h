#ifndef __COMMON
#define __COMMON

enum gender {MALE, FEMALE};

typedef struct _params{
  char *refGenome; /* path to reference genome - fasta */
  char *reps; /* path to repeatmasker file - *rm.out */
  char *dups; /* path to segmental duplications file - bed */
  char *bamFile; /* path to input BAM file */
  char *gaps; /* path to assembly gaps file - bed */
  char *mei;  /* regular expression-like MEI list */
  enum gender sampleGender; /* gender of the sample */
  char runVH; /* boolean stand-in to run VariationHunter */
  char runNS; /* boolean stand-in to run NovelSeq */
  char runSR; /* boolean stand-in to run SPLITREAD */
} parameters;


void init_params(parameters **);
void print_params(parameters *);

#endif
