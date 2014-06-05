#ifndef __COMMANDLINE
#define __COMMANDLINE

#include <getopt.h>
#include <stdlib.h>
#include <stdbool.h>
#include "common.h"

#define FORMAT_MAX_LEN 10
#define INPUT_MAX_LEN 120

#define ERROR_CODE_ARG 10
#define ERROR_CODE_OPTION 11


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
	bool helpWanted;
	bool versionWanted;
} MainOptions;

int parse_command_line( int, char**, parameters*, MainOptions*);
void print_help( void);
void printVersion();
#endif
