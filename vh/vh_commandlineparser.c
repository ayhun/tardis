#include "vh_commandlineparser.h"

#include "stdio.h"
#include "string.h"

char g_error_message[500] = "";

int vh_parseCommand (int argc, char **argv, MainOptions * mainOptions)
{

  int optionIndex = 0;
  char o;
  int errorCode = 0;
  mainOptions->overMapLimit = 1000;	// default value for overlapLimit

  // The first colon at the optstring argument tells getopt to be silent and don't output any text
  while (1)
    {
      static struct option longOptions[] = {
	{"chro", required_argument, NULL, 'c'},
	{"init", required_argument, NULL, 'i'},
	{"prunprob", required_argument, NULL, 'p'},
	{"maxmapping", required_argument, NULL, 'x'},
	{"gap", required_argument, NULL, 'g'},
	{"repeat", required_argument, NULL, 'r'},
	{"svsup", required_argument, NULL, 's'},
	{"format", required_argument, NULL, 'f'},
	{"lib", required_argument, NULL, 'l'},
	{"help", no_argument, NULL, 'h'},
	{"version", no_argument, NULL, 'v'},
	{"output", required_argument, NULL, 'o'},
	{"outputRead", required_argument, NULL, 't'},
	{0, 0, 0, 0}
      };

      o =
	getopt_long_only (argc, argv, ":f:r:o:c:i:g:l:n:t:x:p:s:hv",
			  longOptions, &optionIndex);
      if (o == -1)
	break;

      //TODO: ERROR HANDLING of arguments
      //TODO: TO  BE COMPLETED - Handling input arguments
      switch (o)
	{
	  //TODO: To be completed
	case 'c':
	  strcpy (mainOptions->chroFileName, optarg);
	  break;
	case 'i':
	  strcpy (mainOptions->initializeFileName, optarg);
	  break;
	case 'g':
	  strcpy (mainOptions->gapFileName, optarg);
	  break;
	case 'r':
	  strcpy (mainOptions->repeatFileName, optarg);
	  break;
	case 'o':
	  strcpy (mainOptions->outputFile, optarg);
	  break;
	case 'f':
	  strcpy (mainOptions->format, optarg);
	  break;
	case 'l':
	  strcpy (mainOptions->libFileAdrs, optarg);
	  break;
	case 'p':
	  mainOptions->prunProb = atof (optarg);
	  break;
	case 'x':
	  mainOptions->overMapLimit = atoi (optarg);
	  break;
	case 's':
	  mainOptions->svSup = atof (optarg);
	  break;
	case 't':
	  strcpy (mainOptions->outputRead, optarg);
	  break;
	case 'h':
	  mainOptions->helpWanted = 1;
	  break;
	case 'v':
	  mainOptions->versionWanted = 1;
	  break;
	case ':':		//Missing arguement option
	  errorCode = ERROR_CODE_ARG;
	  sprintf (g_error_message, "Option '%s' requires an argument!",
		   optopt);
	  break;
	  //TODO: Do you want to ignore the invalid options? 
	case '?':		//Invalid option found
	default:
	  errorCode = ERROR_CODE_OPTION;
	  sprintf (g_error_message, "Invalid option found!");
	  break;
	}
    }

  if (optind > argc)
    errorCode = ERROR_CODE_OPTION;

  return errorCode;
}
