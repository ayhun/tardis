#include "stdio.h"
#include "vh_common.h"
#include "vh_repeathandler.h"
#include "vh_divethandler.h"

Repeat *g_repeatTable;
int g_repeatTableSize;
int g_maxRepeatLength = 0;
#define STRMAX 50

int vh_compareInt_Repeat (const void *a, const void *b)
{
  return (((Repeat *) a)->start - ((Repeat *) b)->start);
}


void vh_readRepeatTable (char *repeatFileName)
{

  int satellites;
  char chrom[STRMAX], strand[2], type[STRMAX], class[STRMAX];
  int start, end;
  int res;

  FILE *repeatFile = safe_fopen (repeatFileName, "r");
  if (repeatFile == NULL)
    {
      //TODO: handle me
      fprintf (stderr, "Cannot open file %s\n", repeatFileName);
      exit (-1);
    }

  /* need to implement file format check */
  satellites = 0;
  while (!feof (repeatFile))
    {
      res = fscanf(repeatFile, "%s\t%d\t%d\t%s\t%s\t%s\n", chrom, &start, &end, strand, type, class);

      if (feof (repeatFile))
	break;

      if (res != 6)
	{
	  // todo: global exit function
	  fprintf(stderr, "The repeats file %s is not in the correct BED6 format. Consult the manual on how to generate a valid file.\n", repeatFileName);
	  exit(-1);
	}
    
      if (strstr(class, "Satellite"))
	satellites++;
    }

  if (satellites == 0)
    {
      fprintf(stderr, "No satellites found\n");
    }
  else
    {    
      fprintf(stderr, "Loaded %d satellite regions.\n", satellites);
    }

  g_repeatTable = (Repeat *) malloc(sizeof(Repeat) * satellites);
  
  rewind(repeatFile);

  satellites = 0;
  while (!feof (repeatFile))
    {

      res =  fscanf(repeatFile, "%s\t%d\t%d\t%s\t%s\t%s\n", chrom, &start, &end, strand, type, class);
      if (strstr(class, "Satellite"))
	{
	  g_repeatTable[satellites].chroName = (char *) malloc(sizeof(char) * (strlen(chrom)+1));
	  strcpy(g_repeatTable[satellites].chroName, chrom);
	  g_repeatTable[satellites].start = start;
	  g_repeatTable[satellites].end = end;
	  if (g_maxRepeatLength < g_repeatTable[satellites].end - g_repeatTable[satellites].start)
	    g_maxRepeatLength =  g_repeatTable[satellites].end - g_repeatTable[satellites].start;
	  satellites++;
	}
      
      if (feof (repeatFile))
	break;

      if (res != 6)
	break;

    }
  g_repeatTableSize = satellites;

  qsort (g_repeatTable, g_repeatTableSize, sizeof (Repeat),vh_compareInt_Repeat);


  fclose (repeatFile);
}




int vh_binarySearchInterval (DivetRow * newRow, int pos)	// return 1 if one of the ends falls inside a repeat
{

  int minId = 0, maxId = g_repeatTableSize, midId;
  int stopId, startId;
  int posT = pos - g_maxRepeatLength;
  int i;
  if (posT < 0)
    posT = 0;
  midId = (minId + maxId) / 2;
  while (!
	 (minId > maxId - 3 || g_repeatTable[minId].start == pos
	  || g_repeatTable[maxId].start == pos
	  || g_repeatTable[midId].start == pos))
    {
      midId = (minId + maxId) / 2;
      if (g_repeatTable[midId].start > pos)
	{
	  maxId = midId;

	}
      else if (g_repeatTable[midId].start < pos)
	{
	  minId = midId;
	}
      //      printf("%i %i %i %i %i %i %i\n", minId, maxId, midId, g_repeatTable[minId].start, g_repeatTable[maxId].start, g_repeatTable[midId].start, pos);
    }
  //      printf("%i\n", pos);
  stopId = maxId + 23;
  if (stopId > g_repeatTableSize)
    stopId = g_repeatTableSize;
  startId = 0;
  minId = 0;
  maxId = g_repeatTableSize;
  midId = (minId + maxId) / 2;
  while (!
	 (minId > maxId - 3 || g_repeatTable[minId].end == posT
	  || g_repeatTable[maxId].end == posT
	  || g_repeatTable[midId].end == posT))
    {
      midId = (minId + maxId) / 2;
      if (g_repeatTable[midId].end > posT)
	{
	  maxId = midId;

	}
      else if (g_repeatTable[midId].end < posT)
	{
	  minId = midId;
	}
    }

  startId = minId - 23;
  if (startId < 0)
    startId = 0;
  for (i = startId; i < stopId; i++)
    {
      if (strcmp (g_repeatTable[i].chroName, newRow->chroName) == 0
	  &&
	  ((g_repeatTable[i].end > newRow->locMapLeftEnd
	    && g_repeatTable[i].start < newRow->locMapLeftEnd)
	   || (g_repeatTable[i].end > newRow->locMapRightStart
	       && g_repeatTable[i].start < newRow->locMapRightStart)))
	{
	  return 1;
	}
    }


  //printf("%i %i %i %i %i %i %i\n", startId, stopId, g_repeatTable[startId].start, g_repeatTable[stopId].start, pos, newRow->locMapLeftEnd, newRow->locMapRightStart);
  return 0;




}


int vh_notInRepeat (DivetRow * newRow)
{
  if (vh_binarySearchInterval (newRow, newRow->locMapLeftEnd) || vh_binarySearchInterval (newRow, newRow->locMapRightStart))
    return 0;
  else
    return 1;
}
