#ifndef VH_HASH__H
#define VH_HASH__H

#include<stdlib.h>
#include<string.h>

#define NHASH 29989		//Use a prime number!
#define MULT 31

typedef struct ReadName
{
  char *readName;
  int occurrences;
  double minEdit;
  double sumPhredValue;

  struct ReadName *next;
} ReadName;


//extern ReadName* g_readNameHash[NHASH];

unsigned int vh_hash (char *p);
ReadName *vh_addReadName (ReadName * hash[], char *s, double newEdit,
		       double newPhredValue);
ReadName *getReadNameFromHash (ReadName * hash[], char *readName);
int vh_exportToArray (ReadName * hash[], char *array[], int indexStart);
int vh_countNumReads (ReadName * hash[]);
#endif
