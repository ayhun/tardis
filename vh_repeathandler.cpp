#include "stdio.h"
#include "vh_common.h"
#include "vh_repeathandler.h"
#include "vh_divethandler.h"

Repeat g_repeatTable[500000];
int g_repeatTableSize;
int g_maxRepeatLength=0;

int compareInt_Repeat(const void *a, const void *b)
{
	return (((Repeat*)a)->start - ((Repeat*)b)->start);
}


void readRepeatTable(char* repeatFileName)
{


	FILE* repeatFile = fopen(repeatFileName,"r");
	if (repeatFile == NULL){
		//TODO: handle me
	  fprintf(stderr, "Cannot open file %s\n", repeatFileName);
	  exit(-1);
	}

	int index = 0;
	while (!feof(repeatFile))
	{
		//int res = fscanf(gapFile,"%s%d%d", chroName, &(g_gapTable[index].start), &(g_gapTable[index].end));

		int res = fscanf(repeatFile,"%s%d%d\n", g_repeatTable[index].chroName, &(g_repeatTable[index].start), &(g_repeatTable[index].end));

		if (feof(repeatFile))
		  break;

		//g_gapTable[index].chroName = (char *) malloc((strlen(chroName)+1)*sizeof(char));
		//strcpy(g_gapTable[index].chroName, chroName);
		//printf("RepTable: %s %i %i\n",g_repeatTable[index].chroName, g_repeatTable[index].start, g_repeatTable[index].end);
		
		
		if (res != 3) 
			break;

		if (g_maxRepeatLength<g_repeatTable[index].end-g_repeatTable[index].start)
			g_maxRepeatLength=g_repeatTable[index].end-g_repeatTable[index].start;
		//logDebug(g_gapTable[index].chroName);	
		//g_chroTable[index].chroName = (char *) malloc((strlen(chrName)+1)*sizeof(char));
		//strcpy(g_chroTable[index].chroName, chrName);
	
		index++;
	}
	g_repeatTableSize = index;

	qsort(g_repeatTable, g_repeatTableSize, sizeof(Repeat), compareInt_Repeat);

	
//	for (int i=0; i<g_repeatTableSize; i++)
//		printf("%i\n", g_repeatTable[i].start);
	fclose(repeatFile);
}




bool binarySearchInterval (DivetRow *newRow, int pos) // return true if one of the ends falls inside a repeat
{

int minId=0, maxId=g_repeatTableSize, midId;
int stopId, startId;
int posT=pos-g_maxRepeatLength;
if (posT<0)
	posT=0;
midId=(minId+maxId)/2;
	while (!(minId>maxId-3 || g_repeatTable[minId].start==pos || g_repeatTable[maxId].start==pos || g_repeatTable[midId].start==pos))
	{
		midId=(minId+maxId)/2;
		if (g_repeatTable[midId].start > pos)
		{
			maxId=midId;
			
		} else if (g_repeatTable[midId].start < pos)
		{
			minId=midId;
		}
	//	printf("%i %i %i %i %i %i %i\n", minId, maxId, midId, g_repeatTable[minId].start, g_repeatTable[maxId].start, g_repeatTable[midId].start, pos);
	}
//	printf("%i\n", pos);
stopId=maxId+23;
if(stopId > g_repeatTableSize)
	stopId=g_repeatTableSize;
startId=0;
	minId=0;
	maxId=g_repeatTableSize;
	midId=(minId+maxId)/2;
	while (!(minId>maxId-3 || g_repeatTable[minId].end==posT || g_repeatTable[maxId].end==posT || g_repeatTable[midId].end==posT))
	{
		midId=(minId+maxId)/2;
		if (g_repeatTable[midId].end > posT)
		{
			maxId=midId;
			
		} else if (g_repeatTable[midId].end < posT)
		{
			minId=midId;
		}
	}

startId=minId-23;
if (startId<0)
	startId=0;
for (int i=startId; i<stopId; i++)
{
	if (strcmp(g_repeatTable[i].chroName, newRow->chroName)==0 && ((g_repeatTable[i].end > newRow->locMapLeftEnd && 
		g_repeatTable[i].start < newRow->locMapLeftEnd) || (g_repeatTable[i].end > newRow->locMapRightStart && 
		g_repeatTable[i].start < newRow->locMapRightStart)))
		{
			return true;
		}
}


//printf("%i %i %i %i %i %i %i\n", startId, stopId, g_repeatTable[startId].start, g_repeatTable[stopId].start, pos, newRow->locMapLeftEnd, newRow->locMapRightStart);
return false;




}


bool notInRepeat(DivetRow* newRow)
{
	if (binarySearchInterval(newRow, newRow->locMapLeftEnd)||binarySearchInterval(newRow, newRow->locMapRightStart))
		return false;
	else return true;	
}


