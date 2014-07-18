#ifndef MAIN_H__
#define MAIN_H__

#include "vh_common.h"
#include "vh_commandlineparser.h"
#include "vh_gaphandler.h"

#include <math.h>
#include <string.h>

void vh_printHelp ();
void vh_printVersion ();
//void run(char* fileName, char* chro, char* gap, char *repeat, char* init, double preProsPrune, double minSVSup, char* outputFile);

void vh_clustering (char *, bam_info* , char *, char *, double , char *, char *, int );

#endif
