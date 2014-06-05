#ifndef VH_GAP_HANDLER__H__
#define VH_GAP_HANDLER__H__
#include <stdbool.h>
//Note: Maximum chro name size if considered to be 500
typedef struct Gap {
	char chroName[500];
	//char *chroName;
	int start;
	int end;

} Gap;

extern Gap g_gapTable[1000]; //defined in gaphandler.c	 //TODO: Check the size
extern int g_gapTableSize; //define in gaphandler.c




typedef struct Chro {
	//char *chroName;
	char chroName[500];
	int size;
} Chro;

extern Chro g_chroTable[100]; 
extern int g_chroTableSize;

void readGapTable(char*);
void readInitFile(char*);
void readChros(char*);


#endif

