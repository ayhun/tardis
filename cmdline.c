#include <getopt.h>
#include <stdio.h>
#include <string.h>

#include "tardis.h"
#include "cmdline.h"

int parseCommandLine (int argc, char **argv, parameters *params){

  int index;
  int o;

  static struct option longOptions[] = 
    {
      {"input", required_argument, 0, 'i'},
    };
  
  if (argc == 1){
    printHelp();
    return 0;
  }

  while ( (o = getopt_long ( argc, argv, "i:", longOptions, &index)) != -1 )
    {
      switch (o)
	{
	case 'i':
	  params->bamFile = (char *) malloc ((strlen(optarg)+1) * sizeof(char));
	  strncpy(params->bamFile, optarg, strlen(optarg));
	  break;
	}
    }
}

void printHelp(void){
}

