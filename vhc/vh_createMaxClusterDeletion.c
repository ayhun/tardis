#include "vh_common.h"
#include "vh_hash.h"
#include "vh_divethandler.h"
#include "vh_heap.h"
#include "vh_gaphandler.h"
#include "vh_maximalCluster.h"

//#define NAME_STR_LEN 100


int vh_addToGenomeIndex_Deletion (char *chroName)
{
  LibraryInfo *libInfo;
  DivetRow *divetReadMappingPtr;
  MappingOnGenome *newEl, *newEl2;
  int leftWindowEnd;
  libInfo = g_libInfo;
  while (libInfo != NULL)
    {
      divetReadMappingPtr = libInfo->head;
      while (divetReadMappingPtr != NULL)
	{

	  if (strcmp (divetReadMappingPtr->chroName, chroName) == 0
	      && divetReadMappingPtr->orientationLeft == 'F'
	      && divetReadMappingPtr->orientationRight == 'R'
	      && (divetReadMappingPtr->locMapRightStart -
		  divetReadMappingPtr->locMapLeftEnd > libInfo->maxDelta)
	      && vh_noGap (chroName, divetReadMappingPtr->locMapLeftEnd,
			   divetReadMappingPtr->locMapRightStart)
	      && (divetReadMappingPtr->locMapRightStart -
		  divetReadMappingPtr->locMapLeftEnd < maxDeletionLen))
	    {
	      newEl = (MappingOnGenome *) malloc (sizeof (MappingOnGenome));
	      newEl2 = (MappingOnGenome *) malloc (sizeof (MappingOnGenome));
	      newEl->readMappingPtr = divetReadMappingPtr;
	      newEl2->readMappingPtr = divetReadMappingPtr;
	      newEl->next =
		g_genomeIndexStart[divetReadMappingPtr->locMapLeftEnd];
	      leftWindowEnd =
		divetReadMappingPtr->locMapLeftEnd + libInfo->maxDelta;
	      newEl2->next = g_genomeIndexEnd[leftWindowEnd];
	      g_genomeIndexStart[divetReadMappingPtr->locMapLeftEnd] = newEl;
	      g_genomeIndexEnd[leftWindowEnd] = newEl2;
	    }
	  divetReadMappingPtr = divetReadMappingPtr->next;
	}
      libInfo = libInfo->next;
    }

}

void vh_initializeReadMapping_Deletion (char *chroName, int chroSize)
{
  //Gap info is global    
  LibraryInfo *libInfoPtr = g_libInfo;
  int genomeIndexId;
  int i;

  //Initing the Genome Array - Make sure it is free-ed

  g_genomeIndexStart =
    (MappingOnGenome **) malloc (chroSize * sizeof (MappingOnGenome *));
  g_genomeIndexEnd =
    (MappingOnGenome **) malloc (chroSize * sizeof (MappingOnGenome *));

  if (g_genomeIndexStart == NULL || g_genomeIndexEnd == NULL)
    vh_logWarning ("Memory Problem in vh_createMaxClusterDeletion.cpp:53");
  for (genomeIndexId = 0; genomeIndexId < chroSize; genomeIndexId++)
    {
      g_genomeIndexStart[genomeIndexId] = NULL;
      g_genomeIndexEnd[genomeIndexId] = NULL;
    }

  //Initing the List of begin and end of right side break point ranges - make sure to free
  g_listRightBrkPointIntr =
    (RightBrkPointInterval *) malloc (g_maxListBrkPointIntr *
				      sizeof (RightBrkPointInterval));
  g_tempListRightBrkPointIntr =
    (RightBrkPointInterval *) malloc (g_maxListBrkPointIntr *
				      sizeof (RightBrkPointInterval));
  for (i = 0; i < g_maxListBrkPointIntr; i++)
    {
      g_listRightBrkPointIntr[i].readMappingPtr = NULL;
      //g_listRightBrkPointIntr[i].libInfo=NULL;
      g_tempListRightBrkPointIntr[i].readMappingPtr = NULL;
      //g_tempListRightBrkPointIntr[i].libInfo=NULL;
    }
  g_listRightBrkPointIntrCount = 0;
  vh_addToGenomeIndex_Deletion (chroName);

  /////Malocing the intersectingInterval (intersectinterval) heap
  g_intersectInterval = (Heap *) malloc (sizeof (Heap));
  g_intersectInterval->heapSize = 0;

  while (libInfoPtr != NULL)
    {
      g_maxDeltaAmongLibs = vh_max (g_maxDeltaAmongLibs, libInfoPtr->maxDelta);
      libInfoPtr = libInfoPtr->next;
    }
}


int vh_reevaluate_Deletion (int id, int brkPointLeft)
{
  int brkRightTemp, brkLeftTemp;
  brkRightTemp =
    vh_max (g_tempListRightBrkPointIntr[id].readMappingPtr->locMapRightStart -
	    vh_max (g_tempListRightBrkPointIntr[id].readMappingPtr->
		    libInfo->minDelta - (brkPointLeft -
					 g_tempListRightBrkPointIntr[id].
					 readMappingPtr->locMapLeftEnd), 0),
	    g_tempListRightBrkPointIntr[id].readMappingPtr->locMapLeftEnd);
  brkLeftTemp =
    vh_max (g_tempListRightBrkPointIntr[id].readMappingPtr->locMapRightStart -
	    vh_max (g_tempListRightBrkPointIntr[id].readMappingPtr->
		    libInfo->maxDelta - (brkPointLeft -
					 g_tempListRightBrkPointIntr[id].
					 readMappingPtr->locMapLeftEnd), 0),
	    g_tempListRightBrkPointIntr[id].readMappingPtr->locMapLeftEnd);

  if (g_tempListRightBrkPointIntr[id].keyLorR == 'L')
    {
      g_tempListRightBrkPointIntr[id].locBrkPointLeft = brkLeftTemp;
      g_tempListRightBrkPointIntr[id].key = brkLeftTemp;
      g_tempListRightBrkPointIntr[id].locBrkPointRight = brkRightTemp;
    }
  else if (g_tempListRightBrkPointIntr[id].keyLorR == 'R')
    {
      g_tempListRightBrkPointIntr[id].locBrkPointLeft = brkLeftTemp;
      g_tempListRightBrkPointIntr[id].key = brkRightTemp;
      g_tempListRightBrkPointIntr[id].locBrkPointRight = brkRightTemp;

    }

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int vh_createBreakPointIntervals_Deletion (int brkPointLeft)
{
  MappingOnGenome *ptrMappingOnGenome;
  int newElAdded = 0;
  int genomeId;
  int maxDeltaTemp, minDeltaTemp;
  int locBrkPointLeftTemp, locBrkPointRightTemp;	// These the left and right loci of the interval of the breakpoint on the right
  RightBrkPointInterval *temp;
  int listRightBrkPointIntrId = 0, tempListRightBrkPointIntrId = 0;

  if (g_listRightBrkPointIntrCount > 0)
    {
      ptrMappingOnGenome = g_genomeIndexEnd[brkPointLeft];
      //TODO: Can this be made more efficient using a Heap?
      while (listRightBrkPointIntrId < g_listRightBrkPointIntrCount)
	{
	  if (g_listRightBrkPointIntr
	      [listRightBrkPointIntrId].readMappingPtr->locMapLeftEnd +
	      g_listRightBrkPointIntr
	      [listRightBrkPointIntrId].readMappingPtr->libInfo->maxDelta ==
	      brkPointLeft)
	    {
	      listRightBrkPointIntrId++;
	    }
	  else
	    {
	      vh_copyElBrkPointIntr (tempListRightBrkPointIntrId,listRightBrkPointIntrId);// increaseByOneRightBrkPointIntr(tempListRightBrkPointIntrId);
	      vh_reevaluate_Deletion (tempListRightBrkPointIntrId, brkPointLeft);
	      listRightBrkPointIntrId++;
	      tempListRightBrkPointIntrId++;
	    }
	}
    }


  ptrMappingOnGenome = g_genomeIndexStart[brkPointLeft];
  while (ptrMappingOnGenome != NULL)
    {
      newElAdded = 1;
      locBrkPointRightTemp =
	ptrMappingOnGenome->readMappingPtr->locMapRightStart -
	vh_max ((ptrMappingOnGenome->readMappingPtr->libInfo->minDelta -
		 (brkPointLeft -
		  ptrMappingOnGenome->readMappingPtr->locMapLeftEnd)), 0);
      locBrkPointLeftTemp =
	ptrMappingOnGenome->readMappingPtr->locMapRightStart -
	vh_max ((ptrMappingOnGenome->readMappingPtr->libInfo->maxDelta -
		 (brkPointLeft -
		  ptrMappingOnGenome->readMappingPtr->locMapLeftEnd)), 0);
      g_tempListRightBrkPointIntr
	[tempListRightBrkPointIntrId].locBrkPointRight = locBrkPointRightTemp;
      g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].locBrkPointLeft
	= locBrkPointLeftTemp;
      g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].key =
	locBrkPointLeftTemp;
      g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].readMappingPtr
	= ptrMappingOnGenome->readMappingPtr;
      g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].keyLorR = 'L';
      tempListRightBrkPointIntrId++;
      g_tempListRightBrkPointIntr
	[tempListRightBrkPointIntrId].locBrkPointRight = locBrkPointRightTemp;
      g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].locBrkPointLeft
	= locBrkPointLeftTemp;
      g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].key =
	locBrkPointRightTemp;
      g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].readMappingPtr
	= ptrMappingOnGenome->readMappingPtr;
      g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].keyLorR = 'R';
      tempListRightBrkPointIntrId++;
      ptrMappingOnGenome = ptrMappingOnGenome->next;
    }
  temp = g_listRightBrkPointIntr;
  g_listRightBrkPointIntr = g_tempListRightBrkPointIntr;
  g_tempListRightBrkPointIntr = temp;
  g_listRightBrkPointIntrCount = tempListRightBrkPointIntrId;
  return (newElAdded);
}






// This is the head function called for creating clusters. It startes the loop which iterates the left breakpoint
int vh_createDeletionClusters (int chroSize)
{
  int leftBreakPoint;		// The value of the left breakpoint considered
  // The genome starts from 0
  int startLeftWindow, endLeftWindow;
  g_listRightBrkPointIntrCount = 0;
  int newElAdded = 0;
  // The loop which iterates on the chromosome repersenting the left breakpoint
  for (leftBreakPoint = 1; leftBreakPoint < chroSize; leftBreakPoint++)
    {
      newElAdded = 0;
      newElAdded = vh_createBreakPointIntervals_Deletion (leftBreakPoint);
      if (newElAdded)		// For deletion only when we have added new element we need to check
	vh_createIntersectingIntervals (leftBreakPoint, 2);
    }
  vh_flushOut (g_listPotClusterFound, leftBreakPoint, 2);
}
