#ifndef VH_REPEATHANDLER_H_
#define VH_REPEATHANDLER_H_
#include <stdbool.h>

typedef struct Repeat {
	char chroName[500];
	int start;
	int end;
} Repeat;

extern Repeat g_repeatTable[500000]; //TODO: Malloc instead
extern int g_repeatTableSize;
extern int g_maxRepeatLength;

void readRepeatTable(char*);
bool notInRepeat(DivetRow *);
bool binarySearchInterval(DivetRow *, int);


#endif /*VH_REPEATHANDLER_H_*/
