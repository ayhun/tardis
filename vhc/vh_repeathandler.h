#ifndef VH_REPEATHANDLER_H_
#define VH_REPEATHANDLER_H_

#include "vh_divethandler.h"

typedef struct Repeat
{
  char *chroName;
  int start;
  int end;
} Repeat;

extern struct Repeat *g_repeatTable;	
extern int g_repeatTableSize;
extern int g_maxRepeatLength;

void vh_readRepeatTable (char *);
int vh_notInRepeat (struct DivetRow *);
int vh_binarySearchInterval (struct DivetRow *, int);


#endif /*VH_REPEATHANDLER_H_ */
