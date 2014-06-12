#ifndef _VH_PROCESSING_
#define _VH_PROCESSING_

#include "vh_common.h"
#include "vh_commandlineparser.h"
#include "cmdline.h"
#include "common.h"

#include <math.h>
#include <string.h>
#include <stdbool.h>

void printHelp();	
void printVersion();
//void run(char* fileName, char* chro, char* gap, char *repeat, char* init, double preProsPrune, double minSVSup, char* outputFile);
// parameters* params;
// MainOptions* mainOptions;

#ifdef __cplusplus
extern "C" {
#endif
	void run(char* , char* , char* , char *, char* , double , double , char* ,char*, int);
int vhprocessing(int argc, char** argv, int errorCode);
#ifdef __cplusplus
}
#endif

void pruneAndNormalizeDivets(LibraryInfo* lib, double preProsPrune, int overMapLimit);
void quitProgram(int exitCode);
LibraryInfo* readLibraryInfos(char* libFileAdrs);



// ********************




#endif
