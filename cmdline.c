#include <getopt.h>
#include <stdio.h>
#include <string.h>

#include "tardis.h"
#include "cmdline.h"

int parse_command_line( int argc, char** argv, parameters* params)
{
	int index;
	int o;
	static int is_male = 0, is_female = 0;
	static int run_vh = 0, run_ns = 0, run_sr = 0, run_all = 0;
	static int skip_fastq = 0, skip_sort = 0;

	static struct option long_options[] = 
	{
		{"input"  , required_argument,   0, 'i'},
		{"ref"    , required_argument,   0, 'f'},
		{"gaps"   , required_argument,   0, 'g'},
		{"dups"   , required_argument,   0, 'd'},
		{"reps"   , required_argument,   0, 'r'},      
		{"mei"    , required_argument,   0, 'm'},
		{"threads", required_argument,   0, 't'},
		{"help"   , no_argument,         0, 'h'},
		{"version", no_argument,         0, 'v'},
		{"vh"     , no_argument, &run_vh,     1 },
		{"ns"     , no_argument, &run_ns,     1 },
		{"sr"     , no_argument, &run_sr,     1 },
		{"all"    , no_argument, &run_all,    1 },
		{"xy"     , no_argument,   &is_male,  1 },
		{"xx"     , no_argument, &is_female,  1 },
		{"skip-fastq", no_argument, &skip_fastq,  1 },
		{"skip-sort" , no_argument, &skip_sort,  1 },
		{0        , 0,                   0,  0 }
	};
  
	if( argc == 1)
	{
		print_help();
		return 0;
	}
  
	while( ( o = getopt_long( argc, argv, "hv:i:f:g:d:r:m:", long_options, &index)) != -1)
	  {
		switch( o)
		{
			case 'i':
			  /*
	  			params->bam_file = ( char*) malloc( ( strlen( optarg) + 1) * sizeof( char));
	  			strncpy( params->bam_file, optarg, strlen( optarg));
			  */
	  			set_str( &(params->bam_file), optarg);
	  		break;
	  
			case 'f':
			  /*
				params->ref_genome = ( char*) malloc( ( strlen( optarg) + 1) * sizeof( char));
				strncpy( params->ref_genome, optarg, strlen( optarg));
			  */
				set_str( &(params->ref_genome), optarg);
			break;

			case 'g':
			  /*
				params->gaps = ( char*) malloc( ( strlen( optarg) + 1) * sizeof( char));
				strncpy( params->gaps, optarg, strlen( optarg));
			  */
				set_str( &(params->gaps), optarg);
			break;

			case 'd':
			  /*
				params->dups = ( char*) malloc( ( strlen( optarg) + 1) * sizeof( char));
				strncpy( params->dups, optarg, strlen( optarg));
			  */
			        set_str( &(params->dups), optarg);
			break;

			case 'r':
			  /*
				params->reps = ( char*) malloc( ( strlen( optarg) + 1) * sizeof( char));
				strncpy( params->reps, optarg, strlen( optarg));
			  */
				set_str( &(params->reps), optarg);
			break;

			case 'm':
			  /*
				params->mei = ( char*) malloc( ( strlen( optarg) + 1) * sizeof( char));
				strncpy( params->mei, optarg, strlen( optarg));
			  */
				set_str( &(params->mei), optarg);
			break;

			case 't':
				params->threads = atoi( optarg);
			break;

			case 'h':
				print_help();
				return 0;
			break;

			case 'v':
				fprintf( stdout, "\nTARDIS: Toolkit for Automated and Rapid DIscovery of Structural variants.\n");
				fprintf( stdout, "Version %s. Last update: %s\n\n", VERSION, LAST_UPDATE);
				return 0;
			break; 
		}
	  }
  
	/* TODO: check parameter validity */
  
	/* check algorithms to run; run_all is the default */
	if( !run_vh && !run_sr && !run_ns)
	{
		run_all = 1;
	}
  
	if( run_all)
	{
		run_vh = 1; run_sr = 1; run_ns = 1;
	}
  
	/* check if --xx or --xy is invoked. */
	if( !is_male && !is_female)
	{
		fprintf( stderr, "[TARDIS CMDLINE ERROR] Please select --xx [female] or --xy [male] to specify sample gender.\n");
		return 0;
	}

	if( is_male && is_female)
	{
		fprintf( stderr, "[TARDIS CMDLINE ERROR] Please select either --xx [female] or --xy [male] to specify sample gender. Not both!\n");
		return 0;
	}

	/* check if --input is invoked */
	if( params->bam_file == NULL)
	{
		fprintf( stderr, "[TARDIS CMDLINE ERROR] Please enter input BAM file through the --input option.\n");
		return 0;
	}

	/* check if --ref   is invoked */
	if( params->ref_genome == NULL)
	{
		fprintf( stderr, "[TARDIS CMDLINE ERROR] Please enter reference genome file (FASTA) through the --ref option.\n");
		return 0;
	}

	/* check if --gaps  is invoked */
	if( params->gaps == NULL)
	{
		fprintf( stderr, "[TARDIS CMDLINE ERROR] Please enter the assembly gaps file (BED) through the --gaps option.\n");
		return 0;
	}

	/* check if --reps  is invoked */
	if( params->reps == NULL)
	{
		fprintf( stderr, "[TARDIS CMDLINE ERROR] Please enter the repeats file (RepeaMasker) through the --reps option.\n");
		return 0;
	}

	/* check if --dups  is invoked */
	if( params->dups == NULL)
	{
		fprintf( stderr, "[TARDIS CMDLINE ERROR] Please enter the segmental duplications file (BED) through the --gaps option.\n");
		return 0;
	}

	/* check if --mei   is invoked. If not set Alu:L1Hs:SVA as default */
	if( params->mei == NULL)
	{   
	  /*
   	        params->mei = (char *) malloc(sizeof(char) * (strlen("Alu:L1Hs:SVA")+1));
		strncpy(params->mei, "Alu:L1Hs:SVA", strlen("Alu:L1Hs:SVA"));
	  */
	        set_str( &(params->mei), "Alu:L1Hs:SVA");
	}

	/* check if threads>0 */
	if( params->threads <= 0)
	{
		fprintf( stderr, "[TARDIS CMDLINE WARNING] Invalid number of threads was entered (%d). Resetted to 1.\n", params->threads);
		params->threads = 1;
	}


	/* set flags */
	params->run_vh = run_vh;
	params->run_sr = run_sr;
	params->run_ns = run_ns;
	params->skip_fastq = skip_fastq;
	params->skip_sort = skip_sort;
	if (is_male) params->sample_gender = MALE;
	else if (is_female) params->sample_gender = FEMALE;
}

void print_help( void)
{  
	fprintf( stdout, "\nTARDIS: Toolkit for Automated and Rapid DIscovery of Structural variants.\n");
	fprintf( stdout, "Version %s. Last update: %s\n\n", VERSION, LAST_UPDATE);
	fprintf( stdout, "\t--input [bam file]         : Input file in sorted and indexed BAM format.\n");
	fprintf( stdout, "\t--ref   [reference genome] : Reference genome in FASTA format.\n");
	fprintf( stdout, "\t--gaps  [gaps file]        : Assembly gap coordinates in BED format.\n");
	fprintf( stdout, "\t--dups  [dups file]        : Segmental duplication coordinates in BED format.\n");
	fprintf( stdout, "\t--mei   [\"Alu:L1Hs:SVA\"] : List of mobile element names.\n");
	fprintf( stdout, "\t--xx                       : Sample is male.\n");
	fprintf( stdout, "\t--xy                       : Sample is female.\n");
	fprintf( stdout, "\t--vh                       : Run VariationHunter (read pair + read depth).\n");
	fprintf( stdout, "\t--ns                       : Run NovelSeq (read pair + assembly).\n");
	fprintf( stdout, "\t--sr                       : Run SPLITREAD (split read).\n");
	fprintf( stdout, "\t--all                      : Run all three algorithms above [DEFAULT].\n");
	fprintf( stdout, "\t--skip-fastq               : Skip FASTQ dump for discordants. Use this only if you are regenerating the calls.\n");
	fprintf( stdout, "\t--skip-sort                : Skip FASTQ sort for discordants. Use this only if you are regenerating the calls.\n");
	fprintf( stdout, "\t--version                  : Print version and exit.\n");
	fprintf( stdout, "\t--help                     : Print this help screen and exit.\n\n");
}
