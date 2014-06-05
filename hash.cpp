#ifdef MAIN_HASH
#include <stdio.h>
#endif

#include "vh_hash.h"
#include "vh_logger.h"

//ReadName* g_readNameHash[NHASH];

unsigned int getHash(char *p)
{
	unsigned int h = 0;
	for(; *p; p++)
		h = MULT * h + *p;
	return h % NHASH;
}

ReadName* getReadNameFromHash(ReadName* hash[], char* readName)
{
	int hashedVal = getHash(readName);
	ReadName* p = hash[hashedVal];
	for (; p; p = p->next)
	{
		if (!strcmp(p->readName, readName))
			return p;
	}
	return p;
}

ReadName* addReadName(ReadName* hash[], char *s, double newEdit, double newPhred)
{
	ReadName * p;
	int h = getHash(s);
	for(p = hash[h]; p!= NULL; p = p->next)
	{
		if(strcmp(s, p->readName) == 0)
		{
			p->minEdit = (p->minEdit < newEdit) ? p->minEdit : newEdit;
			p->sumPhredValue += newPhred;
			(p->occurrences)++;
			return p;
		}
	}
	p = (ReadName*)malloc(sizeof(ReadName));
	if(!p){
	//	printf("Out Put Size%i\n", sizeof(ReadName));
		logWarning("Memory Problem");
		//TODO: Memory problem
		return p;
	}
	p->occurrences = 1;
	p->minEdit = newEdit;
	p->sumPhredValue = newPhred;
	p->readName = (char *)malloc((strlen(s)+1)*(sizeof(char)));
	strcpy(p->readName, s);
	p->next = hash[h]; //TODO: Be careful of this line. Shouldn't it be NULL?
	hash[h] = p;

	return p;
}

int exportToArray(ReadName* hash[], char *array[], int indexStart)
{
	ReadName *p;
	for (int count=0; count<NHASH; count++)
	{
		for (p=hash[count]; p!=NULL; p=p->next)
		{
			array[indexStart]=(char*) malloc ((strlen(p->readName)+1)*sizeof(char));
			strcpy(array[indexStart], p->readName);
			indexStart++;
		}
	}

return indexStart;
}

int countNumReads(ReadName *hash[])
{
	int countReads=0;
	ReadName *p;
	for (int count=0; count<NHASH; count++)
	{
		for (p=hash[count]; p!=NULL; p=p->next)
		{
			countReads++;
		}
	}
return countReads;
}
 

#ifdef MAIN_HASH
int main()
{
	char buf[100];
	for (int i = 0; i < NHASH; i++)
		g_readNameHash[i] = NULL;

	float min = 0;
	float mean = 0;
	for (int i = 0; i < 3; i++)
	{
		scanf("%s", buf);
		scanf("%f", &min);
		scanf("%f", &mean);

		addReadName(buf, min, mean);
	}
	for (int i = 0; i < NHASH; i++)
	{
		if (g_readNameHash[i] != NULL)
			printf("\n");
		for (ReadName* p = g_readNameHash[i]; p != NULL; p = p->next)
			printf("|%s %d %f %f|\t", p->readName, p->occurrences, p->minEdit, p->sumPhredValue);
	}	
	printf("\n");
	printf("\n");

	ReadName* p = getReadNameFromHash("sal");
	if(p)
		printf("|%s %d %f %f|\n", p->readName, p->occurrences, p->minEdit, p->sumPhredValue);
	else
		printf("No results found.\n");
	return 0;
}
#endif

