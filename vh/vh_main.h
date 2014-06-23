#ifndef MAIN_H__
#define MAIN_H__

#include "vh_common.h"
#include "vh_commandlineparser.h"

#include <math.h>
#include <string.h>

void vh_printHelp ();
void vh_printVersion ();
//void run(char* fileName, char* chro, char* gap, char *repeat, char* init, double preProsPrune, double minSVSup, char* outputFile);

void vh_run (char *, char *, char *, char *, char *, double, double, char *,
	     char *, int);

#endif
