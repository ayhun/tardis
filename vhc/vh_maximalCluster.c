#include "vh_maximalCluster.h"
#include "vh_common.h"

//RightBrkPointInterval* g_listRightBrkPointIntr;
MappingOnGenome **g_genomeIndexStart;
MappingOnGenome **g_genomeIndexEnd;
RightBrkPointInterval *g_listRightBrkPointIntr;
RightBrkPointInterval *g_tempListRightBrkPointIntr;
int g_maxListBrkPointIntr;
int g_listRightBrkPointIntrCount;
int g_maxDeltaAmongLibs = 0;
ClustersFound *g_listPotClusterFound;
int counter = 0;



int vh_min (int x, int y)
{
  if (x < y)
    return x;
  else
    return y;
}


int vh_max (int x, int y)
{
  if (x < y)
    return y;
  else
    return x;
}

int vh_compare (const void *a, const void *b)
{

  struct RightBrkPointInterval *arg1 = (struct RightBrkPointInterval *) a;
  struct RightBrkPointInterval *arg2 = (struct RightBrkPointInterval *) b;
  if ((*arg1).key > (*arg2).key)
    return 1;
  if ((*arg1).key < (*arg2).key)
    return -1;
  /*	if ((*arg1).key==(*arg1).locBrkPointLeft)
	{
	return 1;
	} else {
	return -1;
	}
  */

  if ((*arg1).key == (*arg2).key)
    {
      if ((*arg1).keyLorR == 'L')
	return -1;
      else
	return 1;
    }
}

int vh_compareReadName (const void *a, const void *b)
{
  return strcmp ((*(DivetRow **) a)->readName->readName,
		 (*(DivetRow **) b)->readName->readName);
}

int vh_compareInt (const void *a, const void *b)
{
  return *(int *) a - *(int *) b;
}

int vh_comparePtr (const void *a, const void *b)
{

  if (strcmp
      ((*(DivetRow **) a)->readName->readName,
       (*(DivetRow **) b)->readName->readName) == 0)
    {
      if ((*(DivetRow **) a)->locMapLeftEnd ==
	  (*(DivetRow **) b)->locMapLeftEnd)
	{
	  return ((*(DivetRow **) a)->locMapRightStart -
		  (*(DivetRow **) b)->locMapRightStart);
	}
      else
	return (*(DivetRow **) a)->locMapLeftEnd -
	  (*(DivetRow **) b)->locMapLeftEnd;
    }
  else
    return strcmp ((*(DivetRow **) a)->readName->readName,
		   (*(DivetRow **) b)->readName->readName);
}



int vh_noGap (char *chroName, int startMapping, int endMapping)
{
  int count;
  for (count = 0; count < g_gapTableSize; count++)
    {
      if (strcmp (chroName, g_gapTable[count].chroName) == 0)
	{
	  if (startMapping < g_gapTable[count].start
	      && endMapping > g_gapTable[count].end)
	    return 0;
	  if (startMapping > g_gapTable[count].start
	      && startMapping < g_gapTable[count].end)
	    return 0;
	  if (endMapping > g_gapTable[count].start
	      && endMapping < g_gapTable[count].end)
	    return 0;
	}
    }
  return 1;
}


void vh_freeLinkedList (MappingOnGenome * cur)
{
  MappingOnGenome *next;
  while (cur != NULL)
    {
      next = cur->next;
      free (cur);
      cur = next;
    }
}

void vh_finalizeReadMapping (char *chroName, int chroSize)
{
  //Free g_genomeIndexStart and g_genomeIndexEnd and their linked list
  int i;
  for (i = 0; i < chroSize; i++)
    {
      vh_freeLinkedList (g_genomeIndexEnd[i]);
      g_genomeIndexEnd[i] = NULL;
      vh_freeLinkedList (g_genomeIndexStart[i]);
      g_genomeIndexStart[i] = NULL;
    }
  free (g_genomeIndexStart);
  g_genomeIndexStart = NULL;
  free (g_genomeIndexEnd);
  g_genomeIndexEnd = NULL;

  //Free g_listRightBrkPointIntrCount and set the g_listRightBrkPointIntrCount
  free (g_listRightBrkPointIntr);
  g_listRightBrkPointIntr = NULL;
  g_listRightBrkPointIntrCount = 0;

  //Free g_tempListRightBrkPointIntr 
  free (g_tempListRightBrkPointIntr);
  g_tempListRightBrkPointIntr = NULL;

  //Free g_intersectInterval -> Heap AND set the heapsize
  free (g_intersectInterval);
  g_intersectInterval = NULL;
}

int vh_copyElBrkPointIntr (int dest, int src)
{
  g_tempListRightBrkPointIntr[dest].locBrkPointLeft =
    g_listRightBrkPointIntr[src].locBrkPointLeft;
  g_tempListRightBrkPointIntr[dest].locBrkPointRight =
    g_listRightBrkPointIntr[src].locBrkPointRight;
  g_tempListRightBrkPointIntr[dest].key = g_listRightBrkPointIntr[src].key;
  g_tempListRightBrkPointIntr[dest].readMappingPtr =
    g_listRightBrkPointIntr[src].readMappingPtr;
  g_tempListRightBrkPointIntr[dest].keyLorR =
    g_listRightBrkPointIntr[src].keyLorR;
}

int vh_isItSubset (int *querySet, int querySetSize, int *patternSet,
		   int patternSetSize)
{
  int idPatternSet = 0;
  int idQuerySet = 0;
  if (patternSetSize == 0)
    return 0;
  while (idQuerySet < querySetSize && idPatternSet < patternSetSize)
    {
      if (querySet[idQuerySet] == patternSet[idPatternSet])
	{
	  idPatternSet++;
	  idQuerySet++;
	}
      else if (querySet[idQuerySet] > patternSet[idPatternSet])
	{
	  idPatternSet++;
	}
      else if (querySet[idQuerySet] < patternSet[idPatternSet])
	return 0;
    }

  if (idPatternSet == patternSetSize && idQuerySet < querySetSize)
    return 0;
  else
    return 1;
}

int vh_flushOut (ClustersFound * listPotClustersFound, int leftBreakPoint,
		 int SVtype)
{
  //printf("%i\n", SVtype);       
  ClustersFound *ptrToOldClusterList;
  ClustersFound *tempPtrOldCluster;
  ptrToOldClusterList = g_listPotClusterFound;
  while (ptrToOldClusterList != NULL)
    {
      if (ptrToOldClusterList->next != NULL
	  && ptrToOldClusterList->next->isMaximalCluster == 0)
	{
	  //              printf("L174\n");
	  tempPtrOldCluster = ptrToOldClusterList->next->next;
	  free (ptrToOldClusterList->next->readMappingPtrArray);
	  free (ptrToOldClusterList->next->readMappingIdArray);
	  free (ptrToOldClusterList->next);
	  ptrToOldClusterList->next = tempPtrOldCluster;
	}
      else if (ptrToOldClusterList->next != NULL
	       && leftBreakPoint >
	       ptrToOldClusterList->next->leftBrkPoint + g_maxDeltaAmongLibs)
	{
	  //              printf("L181\n");
	  vh_outputCluster (ptrToOldClusterList->next, SVtype);
	  tempPtrOldCluster = ptrToOldClusterList->next->next;
	  free (ptrToOldClusterList->next->readMappingPtrArray);
	  free (ptrToOldClusterList->next->readMappingIdArray);
	  free (ptrToOldClusterList->next);
	  ptrToOldClusterList->next = tempPtrOldCluster;
	}
      else
	ptrToOldClusterList = ptrToOldClusterList->next;
    }
  if (g_listPotClusterFound != NULL
      && g_listPotClusterFound->isMaximalCluster == 0)
    {
      tempPtrOldCluster = g_listPotClusterFound->next;
      free (g_listPotClusterFound->readMappingIdArray);
      free (g_listPotClusterFound->readMappingPtrArray);
      free (g_listPotClusterFound);
      g_listPotClusterFound = tempPtrOldCluster;
    }
  else if (g_listPotClusterFound != NULL
	   && leftBreakPoint >
	   g_listPotClusterFound->leftBrkPoint + g_maxDeltaAmongLibs)
    {
      vh_outputCluster (g_listPotClusterFound, SVtype);
      tempPtrOldCluster = g_listPotClusterFound->next;
      free (g_listPotClusterFound->readMappingPtrArray);
      free (g_listPotClusterFound->readMappingIdArray);
      free (g_listPotClusterFound);
      g_listPotClusterFound = tempPtrOldCluster;
    }

}

int vh_addToPotentialOutput (int leftBreakPoint, Heap * heapName, int SVtype)
{
  int count;
  ClustersFound *newCluster;
  ClustersFound *ptrToOldClusterList;
  ClustersFound *tempPtrOldCluster;
  int newClusterIsMaximal = 1;
  int oldClusterIsMaximal = 1;
  newCluster = (ClustersFound *) malloc (sizeof (ClustersFound));
  newCluster->next = NULL;
  newCluster->isMaximalCluster = 1;
  newCluster->leftBrkPoint = leftBreakPoint;
  newCluster->clusterSize = heapName->heapSize;
  newCluster->readMappingPtrArray =
    (DivetRow **) malloc (newCluster->clusterSize * sizeof (DivetRow *));
  newCluster->readMappingIdArray =
    (int *) malloc (newCluster->clusterSize * sizeof (int));
  for (count = 0; count < heapName->heapSize; count++)
    {
      newCluster->readMappingPtrArray[count] =
	heapName->heapArray[count].readMappingPtr;
      newCluster->readMappingIdArray[count] =
	heapName->heapArray[count].readMappingPtr->divetRowId;
    }
  qsort (newCluster->readMappingIdArray, newCluster->clusterSize,sizeof (int), vh_compareInt);
  ptrToOldClusterList = g_listPotClusterFound;
  while (ptrToOldClusterList != NULL && newClusterIsMaximal == 1)
    {
      if (newClusterIsMaximal == 1)
	newClusterIsMaximal =
	  (!vh_isItSubset(newCluster->readMappingIdArray, newCluster->clusterSize,
			  ptrToOldClusterList->readMappingIdArray,
			  ptrToOldClusterList->clusterSize));
      if (newClusterIsMaximal == 0)
	ptrToOldClusterList->leftBrkPoint = leftBreakPoint;
      oldClusterIsMaximal =
	(!vh_isItSubset	 (ptrToOldClusterList->readMappingIdArray,  ptrToOldClusterList->clusterSize, newCluster->readMappingIdArray,  newCluster->clusterSize));
      if (oldClusterIsMaximal == 0 && newClusterIsMaximal == 1)
	ptrToOldClusterList->isMaximalCluster = 0;
      ptrToOldClusterList = ptrToOldClusterList->next;
    }


  if (newClusterIsMaximal == 1)
    {
      newCluster->next = g_listPotClusterFound;
      g_listPotClusterFound = newCluster;
      ptrToOldClusterList = g_listPotClusterFound;
    }
  if (newClusterIsMaximal == 0)
    {
      free (newCluster->readMappingIdArray);
      free (newCluster->readMappingPtrArray);
      free (newCluster);
      newCluster = NULL;
    }

  vh_flushOut (g_listPotClusterFound, leftBreakPoint, SVtype);
}

int vh_outputCluster (ClustersFound * cluster, int SVtype)
{

  int listOfReadsOutputed[10000][3];	// [0]:locMapLeftEnd, [1]: locMapRightStart, [2]: edit distance. This is used to remove the duplicated reads (clonal) as a result of PCR duplication from each cluster.
  int totalAddedToList = 0;
  int clonalRead = 0;
  int readMapCount;
  int countListOutputed;
  //      printf("Cluster size %i\n", cluster->clusterSize);
  if (cluster->clusterSize < 2)
    return 0;
  //              printf("Here\n");
  if (SVtype == 3 && vh_notBothDirections (cluster))
    return 0;
  qsort (cluster->readMappingPtrArray, cluster->clusterSize,
	 sizeof (DivetRow **), vh_compareReadName);
  //      printf("L229\n");
  for (readMapCount = 0; readMapCount < cluster->clusterSize; readMapCount++)
    {
      clonalRead = 0;
      if (readMapCount == 0
	  || strcmp (cluster->readMappingPtrArray[readMapCount]->
		     readName->readName,
		     cluster->readMappingPtrArray[readMapCount -
						  1]->readName->readName) !=
	  0)
	{
	  for (countListOutputed = 0; countListOutputed < totalAddedToList;
	       countListOutputed++)
	    {
	      if (cluster->readMappingPtrArray[readMapCount]->locMapLeftEnd ==
		  listOfReadsOutputed[countListOutputed][0]
		  && cluster->
		  readMappingPtrArray[readMapCount]->locMapRightStart ==
		  listOfReadsOutputed[countListOutputed][1]
		  && cluster->
		  readMappingPtrArray[readMapCount]->editDistance ==
		  listOfReadsOutputed[countListOutputed][2])
		clonalRead = 1;
	    }
	  if (clonalRead == 0)
	    {
	      fprintf (fileOutput, "%s %s %i %s %i %i %g %g %s %c %c ",
		       cluster->readMappingPtrArray[readMapCount]->
		       readName->readName,
		       cluster->readMappingPtrArray[readMapCount]->chroName,
		       cluster->
		       readMappingPtrArray[readMapCount]->locMapLeftEnd,
		       cluster->readMappingPtrArray[readMapCount]->chroName,
		       cluster->
		       readMappingPtrArray[readMapCount]->locMapRightStart,
		       SVtype,
		       cluster->readMappingPtrArray[readMapCount]->phredScore,
		       cluster->
		       readMappingPtrArray[readMapCount]->editDistance,
		       cluster->readMappingPtrArray[readMapCount]->
		       libInfo->libName,
		       cluster->
		       readMappingPtrArray[readMapCount]->orientationLeft,
		       cluster->
		       readMappingPtrArray[readMapCount]->orientationRight);
	      listOfReadsOutputed[totalAddedToList][0] =
		cluster->readMappingPtrArray[readMapCount]->locMapLeftEnd;
	      listOfReadsOutputed[totalAddedToList][1] =
		cluster->readMappingPtrArray[readMapCount]->locMapRightStart;
	      listOfReadsOutputed[totalAddedToList][2] =
		(int) cluster->
		readMappingPtrArray[readMapCount]->editDistance;
	      totalAddedToList++;

	    }
	}
    }
  //printf("Done\n");
  fprintf (fileOutput, "END\n");
}

int vh_createIntersectingIntervals (int leftBreakPoint, int SVtype)
{
  //      printf("l:269\n");
  int newElAdded, countIntrEndPoints;
  newElAdded = 0;
  qsort (g_listRightBrkPointIntr, g_listRightBrkPointIntrCount,
	 sizeof (struct RightBrkPointInterval), vh_compare);
  for (countIntrEndPoints = 0;
       countIntrEndPoints < g_listRightBrkPointIntrCount;
       countIntrEndPoints++)
    {
      counter++;
      //      if (counter%10000==0)
      //              printf("L302 %i %i %i\n", g_listRightBrkPointIntrCount, g_intersectInterval->heapSize, leftBreakPoint);

      //if (g_listRightBrkPointIntr[countIntrEndPoints].locBrkPointLeft==g_listRightBrkPointIntr[countIntrEndPoints].key)
      if (g_listRightBrkPointIntr[countIntrEndPoints].keyLorR == 'L')
	{
	  newElAdded = 1;
	  vh_addToHeap (g_listRightBrkPointIntr
			[countIntrEndPoints].readMappingPtr,
			g_listRightBrkPointIntr
			[countIntrEndPoints].locBrkPointRight,
			g_intersectInterval);

	  //if (leftBreakPoint>500 && leftBreakPoint<1500)
	  //      printf("%i %i\n", leftBreakPoint, g_intersectInterval->heapSize);

	  //for (int testCount=0; testCount<g_intersectInterval->heapSize; testCount++)
	  //{
	  //}
	  //writeHeap(g_intersectInterval);

	}
      else if (g_listRightBrkPointIntr[countIntrEndPoints].keyLorR == 'R')
	//(g_listRightBrkPointIntr[countIntrEndPoints].locBrkPointRight==g_listRightBrkPointIntr[countIntrEndPoints].key)
	{
	  if (vh_minValue_heap (g_intersectInterval) ==
	      g_listRightBrkPointIntr[countIntrEndPoints].key)
	    {
	      if (newElAdded == 1)
		{
		  vh_addToPotentialOutput (leftBreakPoint, g_intersectInterval,
					   SVtype);
		  newElAdded = 0;
		}
	      //writeHeap(g_intersectInterval);
	      vh_heap_remove_top (g_intersectInterval);
	    }
	  else
	    {
	      printf ("An Error Occured\n");
	      vh_writeHeap (g_intersectInterval);
	    }
	}
    }
  /*if (leftBreakPoint>500 && leftBreakPoint<1500)
    printf("%i %i\n", leftBreakPoint, g_intersectInterval->heapSize);
  */
  g_intersectInterval->heapSize = 0;
}

int vh_notBothDirections (ClustersFound * cluster)
{
  int FF, RR, readMapCount;
  FF = 0;
  RR = 0;
  for (readMapCount = 0; readMapCount < cluster->clusterSize; readMapCount++)
    {
      if (cluster->readMappingPtrArray[readMapCount]->orientationLeft == 'F'
	  && cluster->readMappingPtrArray[readMapCount]->orientationRight ==
	  'F')
	FF = 1;

      if (cluster->readMappingPtrArray[readMapCount]->orientationLeft == 'R'
	  && cluster->readMappingPtrArray[readMapCount]->orientationRight ==
	  'R')
	RR = 1;

      if (FF && RR)
	return 0;
    }
  return 1;
}
