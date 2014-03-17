#include <getopt.h>
#include <stdio.h>
#include <string.h>

#include "tardis.h"
#include "cmdline.h"

int parseCommandLine (int argc, char **argv, parameters *params){

  int index;
  int o;
  static int isMale=0, isFemale=0;
  static int runVH=0, runNS=0, runSR=0, runAll=0;
  
  static struct option longOptions[] = 
    {
      {"input"  , required_argument,   0, 'i'},
      {"ref"    , required_argument,   0, 'f'},
      {"gaps"   , required_argument,   0, 'g'},
      {"dups"   , required_argument,   0, 'd'},
      {"reps"   , required_argument,   0, 'r'},      
      {"mei"    , required_argument,   0, 'm'},
      {"help"   , no_argument,         0, 'h'},
      {"version", no_argument,         0, 'v'},
      {"vh"     , no_argument, &runVH,     1 },
      {"ns"     , no_argument, &runNS,     1 },
      {"sr"     , no_argument, &runSR,     1 },
      {"all"    , no_argument, &runAll,    1 },
      {"xy"     , no_argument,   &isMale,  1 },
      {"xx"     , no_argument, &isFemale,  1 },
      {0        , 0,                   0,  0 }
    };
  
  if (argc == 1){
    printHelp();
    return 0;
  }
  
  while ( (o = getopt_long ( argc, argv, "hv:i:f:g:d:r:m:", longOptions, &index)) != -1 )
    {
      switch (o)
	{
	case 'i':
	  params->bamFile = (char *) malloc ((strlen(optarg)+1) * sizeof(char));
	  strncpy(params->bamFile, optarg, strlen(optarg));
	  break;
	  
	case 'f':
	  params->refGenome = (char *) malloc ((strlen(optarg)+1) * sizeof(char));
	  strncpy(params->refGenome, optarg, strlen(optarg));
	  break;
	  
	case 'g':
	  params->gaps = (char *) malloc ((strlen(optarg)+1) * sizeof(char));
	  strncpy(params->gaps, optarg, strlen(optarg));
	  break;
	  
	case 'd':
	  params->dups = (char *) malloc ((strlen(optarg)+1) * sizeof(char));
	  strncpy(params->dups, optarg, strlen(optarg));
	  break;
	  
	case 'r':
	  params->reps = (char *) malloc ((strlen(optarg)+1) * sizeof(char));
	  strncpy(params->reps, optarg, strlen(optarg));
	  break;
	  
	case 'm':
	  params->mei = (char *) malloc ((strlen(optarg)+1) * sizeof(char));
	  strncpy(params->mei, optarg, strlen(optarg));
	  break;
	  
	case 'h':
	  printHelp();
	  return 0;
	  break;
	  
	case 'v':
	  fprintf(stdout, "\nTARDIS: Toolkit for Automated and Rapid DIscovery of Structural variants.\n");
	  fprintf(stdout, "Version %s. Last update: %s\n\n", VERSION, LAST_UPDATE);
	  return 0;
	  break;
	  
	}
    }
  
  /* check parameter validity */
  
  /* check algorithms to run */
  if (!runVH && !runSR && !runNS) // runAll is the default
    runAll = 1;
  
  if (runAll){
    runVH = 1; runSR = 1; runNS = 1;
  }
  
  /* check if --xx or --xy is invoked. */
  
  if (!isMale && !isFemale){
    fprintf(stderr, "[TARDIS CMDLINE ERROR] Please select --xx [female] or --xy [male] to specify sample gender.\n");
    return 0;
  }

  if (isMale && isFemale){
    fprintf(stderr, "[TARDIS CMDLINE ERROR] Please select either --xx [female] or --xy [male] to specify sample gender. Not both!\n");
    return 0;
  }

  /* check if --input is invoked */
  if (params->bamFile == NULL){
    return 0;
  }

  /* check if --ref   is invoked */
  if (params->refGenome == NULL){
    return 0;
  }

  /* check if --gaps  is invoked */
  if (params->gaps == NULL){
    return 0;
  }

  /* check if --reps  is invoked */
  if (params->reps == NULL){
    return 0;
  }

  /* check if --dups  is invoked */
  if (params->dups == NULL){
    return 0;
  }

  /* check if --mei   is invoked -- this can be optional */
  
  if (params->mei == NULL){
    return 0;
  }


}

void printHelp(void){
  
  fprintf(stdout, "\nTARDIS: Toolkit for Automated and Rapid DIscovery of Structural variants.\n");
  fprintf(stdout, "Version %s. Last update: %s\n\n", VERSION, LAST_UPDATE);
  fprintf(stdout, "\t--input [bam file]        : Input file in sorted and indexed BAM format.\n");
  fprintf(stdout, "\t--ref   [reference genome]: Reference genome in FASTA format.\n");
  fprintf(stdout, "\t--gaps  [gaps file]       : Assembly gap coordinates in BED format.\n");
  fprintf(stdout, "\t--dups  [dups file]       : Segmental duplication coordinates in BED format.\n");
  fprintf(stdout, "\t--mei   [\"Alu|L1\"]      : List of mobile element names.\n");
  fprintf(stdout, "\t--xx                      : Sample is male.\n");
  fprintf(stdout, "\t--xy                      : Sample is female.\n");
  fprintf(stdout, "\t--vh                      : Run VariationHunter (read pair + read depth).\n");
  fprintf(stdout, "\t--ns                      : Run NovelSeq (read pair + assembly).\n");
  fprintf(stdout, "\t--sr                      : Run SPLITREAD (split read).\n");
  fprintf(stdout, "\t--all                     : Run all three algorithms above [DEFAULT].\n");
  fprintf(stdout, "\t--version                 : Print version and exit.\n");
  fprintf(stdout, "\t--help                    : Print this help screen and exit.\n\n");
}

