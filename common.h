#ifndef __COMMON
#define __COMMON

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <zlib.h>

/* Exit Codes */
#define EXIT_SUCCESS 0
#define EXIT_COMMON 1
#define EXIT_MAXBAMS 2
#define EXIT_PARAM_ERROR 3
#define EXIT_EXTERNAL_PROG_ERROR 4
#define EXIT_FILE_OPEN_ERROR 5

/* Return Codes */
#define RETURN_SUCCESS 0
#define RETURN_ERROR 1

#define MAX_BAMS 256

// Track memory usage
extern long long memUsage;

enum gender{ MALE, FEMALE};

typedef struct _params
{
	char* ref_genome; /* path to reference genome - fasta */
	char* reps; /* path to repeatmasker file - *rm.out */
	char* dups; /* path to segmental duplications file - bed */
  	char* bam_files; /* paths to comma separated input BAM files as a single string before being tokenized */
	char* bam_list_path; /* path to a file that lists BAM file paths in advance */
	char** bam_file_list; /* the actual list that holds all bam file paths after tokenization */
	char* gaps; /* path to assembly gaps file - bed */
	char* mei;  /* regular expression-like MEI list */
        char* outprefix; /* prefix for the output files */
  /* move the gender to a new samples data structure, and add a pointer to sample 
        enum gender sample_gender; /* gender of the sample 
  */
        int  force_read_length; /* force read length to a certain value, discard those that are shorter. Hidden feature due to GIAB */
	char run_vh; /* boolean stand-in to run VariationHunter */
	char run_ns; /* boolean stand-in to run NovelSeq */
	char run_sr; /* boolean stand-in to run SPLITREAD */
	char skip_fastq; /* boolean stand-in to skip FASTQ dump */
	char skip_sort; /* boolean stand-in to skip FASTQ sort */
	char skip_remap; /* boolean stand-in to skip FASTQ remap */
	char skip_vhcluster; /* boolean stand-in to skip VH clustering */
	int  threads; /* number of threads to use for parallel mrFAST, and maybe future parallelization of TARDIS */
	int num_bams; /* number of input BAM files */
} parameters;

/* Parameter related TARDIS functions */
void init_params( parameters**);
void print_params( parameters*);

/* FILE opening and error printing functions. For opening regular and BAM/SAM
 files safely */
void print_error( char*);
FILE* safe_fopen( char* path, char* mode);
gzFile safe_fopen_gz( char* path, char* mode);
htsFile* safe_hts_open( char* path, char* mode);

/* General BAM processing functions */
int is_proper( int flag);
int is_concordant( bam1_core_t bam_alignment_core, int min, int max);
char base_as_char( int base_as_int);
char complement_char( char base);
void qual_to_ascii( char* qual);

/* String functions */
void set_str( char **target, char *source); /* Even safer than strncpy */
void reverse_string( char* str);

/* Misc. Utility */
int compare_size_int( const void* p, const void* q);
void print_quote( void);

// Memory allocation/tracking functions
void* getMem( size_t size);
void freeMem( void* ptr, size_t size);
double getMemUsage();

#endif
