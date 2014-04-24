#ifndef __PROCESSBAM
#define __PROCESSBAM

/* htslib headers */
#include <htslib/sam.h>
#include <htslib/hts.h>

typedef struct _bam_info
{
	htsFile* bam_file; /* file pointer to the BAM file */
  	int num_chrom; /* number of chromosomes */
	int* chrom_lengths; /* lengths of the chromosomes */
	char** chrom_names; /* names of the chromosomes */
  	char* sample_name; /* name of the sample, parsed from SM in the BAM header */
	float frag_avg; /* average fragment size */
  	float frag_std; /* fragment size standard deviation */
	int frag_med; /* median of the fragment sizes */
  	int conc_min; /* min cutoff for concordants */
  	int conc_max; /* max cutoff for concordants */
} bam_info;

/* Function Prototypes */
void load_bam( bam_info* in_bam, char* path);
void print_bam( bam_info* in_bam);
char base_as_char( int base_as_int);
void get_sample_name( bam_info* in_bam, char* header_text);
int compare_size( const void* p, const void* q);

#endif
