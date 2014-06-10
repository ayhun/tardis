#ifndef VH_HASH__H
#define VH_HASH__H

#include<stdlib.h>
#include<string.h>
#include <stdbool.h>

#define NHASH 29989 //Use a prime number!
#define MULT 31

struct ReadNameS
{
	char* readName;
	int occurrences;
	double minEdit;
	double sumPhredValue;

	struct ReadNameS* next;
};

typedef struct ReadNameS ReadName;

//extern ReadName* g_readNameHash[NHASH];

unsigned int hash(char *p);
ReadName* addReadName(ReadName* hash[], char *s, double newEdit, double newPhredValue);
ReadName* getReadNameFromHash(ReadName* hash[], char* readName);
int exportToArray(ReadName* hash[], char *array[], int indexStart);
int countNumReads(ReadName *hash[]);
#endif

