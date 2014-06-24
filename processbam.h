#ifndef __PROCESSBAM
#define __PROCESSBAM

/* htslib headers */
#include <htslib/sam.h>
#include <htslib/hts.h>
#include "common.h"

struct library_properties
{
	char* libname; /* id/name of the library */
	float frag_avg; /* average fragment size */
  	float frag_std; /* fragment size standard deviation */
	int frag_med; /* median of the fragment sizes */
  	int conc_min; /* min cutoff for concordants */
  	int conc_max; /* max cutoff for concordants */
	char* fastq1; /* file name for the FASTQ file of the /1 reads */
	char* fastq2; /* file name for the FASTQ file of the /2 reads */
};

typedef struct _bam_info
{
	htsFile* bam_file; /* file pointer to the BAM file */
  	int num_chrom; /* number of chromosomes */
	int* chrom_lengths; /* lengths of the chromosomes */
	char** chrom_names; /* names of the chromosomes */
  	char* sample_name; /* name of the sample, parsed from SM in the BAM header */
	int num_libraries; /* number of libraries, counted from the RG tags in the BAM header */
	struct library_properties** libraries; /* each library_properties struct holds statistical/other info */
} bam_info;

/* Function Prototypes */
void load_bam( bam_info* in_bam, char* path);
void print_bam( bam_info* in_bam);
void print_libs( bam_info* in_bam);
int find_library_index( bam_info* in_bam, char* library_name);
int sufficient_fragments_sampled( int* fragments_sampled, int num_libraries);
void create_fastq( bam_info* in_bam, parameters *params);

/* BAM Utility functions */
void get_sample_name( bam_info* in_bam, char* header_text);
void get_library_count( bam_info* in_bam, char* header_text);
void get_library_names( bam_info* in_bam, char* header_text);

#endif
