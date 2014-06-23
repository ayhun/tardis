#ifndef __COMMON
#define __COMMON

#include <htslib/sam.h>
#define FORMAT_MAX_LEN 10
#define INPUT_MAX_LEN 120
#define ERROR_CODE_ARG 10
#define ERROR_CODE_OPTION 11



enum gender{ MALE, FEMALE};

typedef struct _params
{
	char* ref_genome; /* path to reference genome - fasta */
	char* reps; /* path to repeatmasker file - *rm.out */
	char* dups; /* path to segmental duplications file - bed */
	char* bam_file; /* path to input BAM file */
	char* gaps; /* path to assembly gaps file - bed */
	char* mei;  /* regular expression-like MEI list */
	enum gender sample_gender; /* gender of the sample */
	char run_vh; /* boolean stand-in to run VariationHunter */
	char run_ns; /* boolean stand-in to run NovelSeq */
	char run_sr; /* boolean stand-in to run SPLITREAD */
	char skip_fastq; /* boolean stand-in to skip FASTQ dump */
	char skip_sort; /* boolean stand-in to skip FASTQ sort */
	int  threads; /* number of threads to use for parallel mrFAST, and maybe future parallelization of TARDIS */
} parameters;

typedef struct MainOptions
{
	float prunProb;
	double svSup;
	char format[FORMAT_MAX_LEN]; //TODO: Can be removed
	char libFileAdrs[INPUT_MAX_LEN];
	char chroFileName[INPUT_MAX_LEN];
	char repeatFileName[INPUT_MAX_LEN];
	char initializeFileName[INPUT_MAX_LEN];
	char gapFileName[INPUT_MAX_LEN];
	char outputFile[INPUT_MAX_LEN];
	char outputRead[INPUT_MAX_LEN];
	int overMapLimit;
	int helpWanted;
	int versionWanted;
} MainOptions;

/* Function prototypes */
void init_params( parameters**);
void print_params( parameters*);
void print_error( char*);
FILE* gfOpen( char*, char*);
int is_concordant( bam1_core_t bam_alignment_core, int min, int max);
char complement_char( char base);
void reverse_string( char* str);
void set_str( char **target, char *source);

#endif
