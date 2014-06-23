#ifndef TEST_GET_OPT_H__
#define TEST_GET_OPT_H__

#include <getopt.h>
#include <stdlib.h>

#define FORMAT_MAX_LEN 10
#define INPUT_MAX_LEN 120

#define ERROR_CODE_ARG 10
#define ERROR_CODE_OPTION 11

typedef struct MainOptions
{
  float prunProb;
  double svSup;
  char format[FORMAT_MAX_LEN];	//TODO: Can be removed
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

int vh_parseCommand (int argc, char **argv, MainOptions * mainOptions);
void vh_printHelp ();
void vh_printVersion ();

#endif
