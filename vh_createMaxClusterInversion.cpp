#include "vh_common.h"
#include "vh_hash.h"
#include "vh_divethandler.h"
#include "vh_heap.h"
#include "vh_gaphandler.h"
#include "vh_maximalCluster.h"

//#define NAME_STR_LEN 100



int addToGenomeIndex_Inversion(char *chroName)
{
	LibraryInfo* libInfo;
	DivetRow* divetReadMappingPtr;
	MappingOnGenome *newEl, *newEl2;
	int leftWindowEnd, leftWindowStart;	
	libInfo=g_libInfo;
	while (libInfo!=NULL)
	{
		divetReadMappingPtr=libInfo->head;
		while(divetReadMappingPtr!=NULL)
		{
			
			if (strcmp(divetReadMappingPtr->chroName, chroName)==0 && ((divetReadMappingPtr->orientationLeft=='F' && divetReadMappingPtr->orientationRight=='F') || (divetReadMappingPtr->orientationLeft=='R' && divetReadMappingPtr->orientationRight=='R' )) && (divetReadMappingPtr->locMapLeftEnd<divetReadMappingPtr->locMapRightStart) && (divetReadMappingPtr->locMapRightStart - divetReadMappingPtr->locMapLeftEnd<maxInversionLen )) 
			{


				newEl=(MappingOnGenome *) malloc(sizeof(MappingOnGenome));
				newEl2=(MappingOnGenome *) malloc(sizeof(MappingOnGenome));
				newEl->readMappingPtr=divetReadMappingPtr;
				newEl2->readMappingPtr=divetReadMappingPtr;
				if (divetReadMappingPtr->orientationLeft=='F' && divetReadMappingPtr->orientationRight=='F')
				{
					newEl->next=g_genomeIndexStart[divetReadMappingPtr->locMapLeftEnd];
					leftWindowEnd=divetReadMappingPtr->locMapLeftEnd+libInfo->maxDelta;
					newEl2->next=g_genomeIndexEnd[leftWindowEnd];
					g_genomeIndexStart[divetReadMappingPtr->locMapLeftEnd]=newEl;
					g_genomeIndexEnd[leftWindowEnd]=newEl2;

				}
				else if (divetReadMappingPtr->orientationLeft=='R' && divetReadMappingPtr->orientationRight=='R')
				{
					leftWindowStart=max(divetReadMappingPtr->locMapLeftStart-libInfo->maxDelta,0);
					newEl->next=g_genomeIndexStart[leftWindowStart];
					newEl2->next=g_genomeIndexEnd[divetReadMappingPtr->locMapLeftStart];
					g_genomeIndexStart[leftWindowStart]=newEl;
					g_genomeIndexEnd[divetReadMappingPtr->locMapLeftStart]=newEl2;

				}
				/*
				newEl->next=g_genomeIndexStart[divetReadMappingPtr->locMapLeftEnd];
				leftWindowEnd=divetReadMappingPtr->locMapLeftEnd+libInfo->maxDelta;
				newEl2->next=g_genomeIndexEnd[leftWindowEnd];
				g_genomeIndexStart[divetReadMappingPtr->locMapLeftEnd]=newEl;
				g_genomeIndexEnd[leftWindowEnd]=newEl2;
				*/
			}
			divetReadMappingPtr=divetReadMappingPtr->next;
		}
		libInfo=libInfo->next;
	}

}

void initializeReadMapping_Inversion(char* chroName, int chroSize)
{
	//Gap info is global	
	LibraryInfo* libInfoPtr=g_libInfo;
	//Initing the Genome Array - Make sure it is free-ed
	
	g_genomeIndexStart = (MappingOnGenome**) malloc (chroSize *sizeof(MappingOnGenome*));
	g_genomeIndexEnd = (MappingOnGenome**) malloc (chroSize * sizeof(MappingOnGenome*));

	if (g_genomeIndexStart == NULL || g_genomeIndexEnd == NULL)
		logWarning("Memory Problem");
	for (int genomeIndexId = 0; genomeIndexId < chroSize; genomeIndexId++)
	{
		g_genomeIndexStart[genomeIndexId] = NULL;
		g_genomeIndexEnd[genomeIndexId] = NULL;
	}
	printf("Inversion\n");
	//Initing the List of begin and end of right side break point ranges - make sure to free
	g_listRightBrkPointIntr = (RightBrkPointInterval*) malloc (g_maxListBrkPointIntr * sizeof(RightBrkPointInterval));
	g_tempListRightBrkPointIntr=(RightBrkPointInterval*) malloc (g_maxListBrkPointIntr * sizeof(RightBrkPointInterval));
	for (int i=0; i<g_maxListBrkPointIntr; i++)
	{
		g_listRightBrkPointIntr[i].readMappingPtr=NULL;
		//g_listRightBrkPointIntr[i].libInfo=NULL;
		g_tempListRightBrkPointIntr[i].readMappingPtr=NULL;
		//g_tempListRightBrkPointIntr[i].libInfo=NULL;
	}
	g_listRightBrkPointIntrCount = 0;
	addToGenomeIndex_Inversion(chroName);

	/////Malocing the intersectingInterval (intersectinterval) heap
	g_intersectInterval = (Heap*) malloc(sizeof(Heap));
	g_intersectInterval->heapSize = 0;
	
	while(libInfoPtr!=NULL)
	{
		g_maxDeltaAmongLibs=max(g_maxDeltaAmongLibs, libInfoPtr->maxDelta);
		libInfoPtr=libInfoPtr->next;
	}
}
///////////////////////////////////The left end side and right end side of the Right breakpoint//////////////
int reevaluate_rightBrkPoint_FF(int id, int brkPointLeft)
{
int brkRightTemp, brkLeftTemp;
// we need to take the max with locMapRightStart, because it might be possible that (brkPointLeft-g_tempListRightPointIntr[id].readMappingPtr->locMapLeftEnd) be bigger than minDelta.
brkRightTemp = max (g_tempListRightBrkPointIntr[id].readMappingPtr->locMapRightEnd + 
		    g_tempListRightBrkPointIntr[id].readMappingPtr->libInfo->maxDelta-
		   (brkPointLeft - g_tempListRightBrkPointIntr[id].readMappingPtr->locMapLeftEnd)
	  	   , g_tempListRightBrkPointIntr[id].readMappingPtr->locMapRightEnd);

brkLeftTemp = max ( g_tempListRightBrkPointIntr[id].readMappingPtr->locMapRightEnd + 
		    g_tempListRightBrkPointIntr[id].readMappingPtr->libInfo->minDelta-
		    (brkPointLeft - g_tempListRightBrkPointIntr[id].readMappingPtr->locMapLeftEnd), 
		    g_tempListRightBrkPointIntr[id].readMappingPtr->locMapRightEnd);

if (g_tempListRightBrkPointIntr[id].keyLorR=='L')
{
	g_tempListRightBrkPointIntr[id].locBrkPointLeft = brkLeftTemp;
	g_tempListRightBrkPointIntr[id].key=brkLeftTemp;
	g_tempListRightBrkPointIntr[id].locBrkPointRight=brkRightTemp;
} else if (g_tempListRightBrkPointIntr[id].keyLorR=='R')
{
	g_tempListRightBrkPointIntr[id].locBrkPointLeft = brkLeftTemp;
	g_tempListRightBrkPointIntr[id].key=brkRightTemp;
	g_tempListRightBrkPointIntr[id].locBrkPointRight=brkRightTemp;
	
}


}

int reevaluate_rightBrkPoint_RR(int id, int brkPointLeft)
{
int brkRightTemp, brkLeftTemp;
brkRightTemp = min(max(g_tempListRightBrkPointIntr[id].readMappingPtr->locMapRightStart -
		  (g_tempListRightBrkPointIntr[id].readMappingPtr->libInfo->minDelta - 
		  (g_tempListRightBrkPointIntr[id].readMappingPtr->locMapLeftEnd - brkPointLeft)),
		   g_tempListRightBrkPointIntr[id].readMappingPtr->locMapLeftStart), g_tempListRightBrkPointIntr[id].readMappingPtr->locMapRightStart);

brkLeftTemp = max(g_tempListRightBrkPointIntr[id].readMappingPtr->locMapRightStart -
		 (g_tempListRightBrkPointIntr[id].readMappingPtr->libInfo->maxDelta - 
		 (g_tempListRightBrkPointIntr[id].readMappingPtr->locMapLeftStart - brkPointLeft)),
		  g_tempListRightBrkPointIntr[id].readMappingPtr->locMapLeftStart);

if (g_tempListRightBrkPointIntr[id].keyLorR=='L')
{
	g_tempListRightBrkPointIntr[id].locBrkPointLeft = brkLeftTemp;
	g_tempListRightBrkPointIntr[id].key=brkLeftTemp;
	g_tempListRightBrkPointIntr[id].locBrkPointRight=brkRightTemp;
} else if (g_tempListRightBrkPointIntr[id].keyLorR=='R')
{
	g_tempListRightBrkPointIntr[id].locBrkPointLeft = brkLeftTemp;
	g_tempListRightBrkPointIntr[id].key=brkRightTemp;
	g_tempListRightBrkPointIntr[id].locBrkPointRight=brkRightTemp;
	
}


}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool createBreakPointIntervals_Inversion(int brkPointLeft)
{
	MappingOnGenome* ptrMappingOnGenome;
	bool newElAdded = false;
	int genomeId;
	int maxDeltaTemp, minDeltaTemp;
	int locBrkPointLeftTemp, locBrkPointRightTemp; // These the left and right loci of the interval of the breakpoint on the right
	RightBrkPointInterval *temp;
	int listRightBrkPointIntrId = 0, tempListRightBrkPointIntrId = 0;

	if (g_listRightBrkPointIntrCount > 0)
	{
		ptrMappingOnGenome = g_genomeIndexEnd[brkPointLeft];
		//TODO: Can this be made more efficient using a Heap?
		while (listRightBrkPointIntrId < g_listRightBrkPointIntrCount)
		{

			if (g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->orientationLeft=='F' && g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->orientationRight=='F')
			{
				if ((g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->locMapLeftEnd + g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->libInfo->maxDelta == brkPointLeft) || (brkPointLeft==g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->locMapRightEnd))
				{	
					listRightBrkPointIntrId++;
				}
				else
				{
					copyElBrkPointIntr(tempListRightBrkPointIntrId, listRightBrkPointIntrId);
					//increaseByOneRightBrkPointIntr(tempListRightBrkPointIntrId);
					//decreaseByOneRightBrkPointIntr_FF(tempListRightBrkPointIntrId);
					reevaluate_rightBrkPoint_FF(tempListRightBrkPointIntrId, brkPointLeft);
					listRightBrkPointIntrId++;
					tempListRightBrkPointIntrId++;
				}
			} else if (g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->orientationLeft=='R' && g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->orientationRight=='R')
			{
				if ((g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->locMapLeftStart ==  brkPointLeft) || (g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->locMapRightStart-brkPointLeft < g_listRightBrkPointIntr[listRightBrkPointIntrId].readMappingPtr->libInfo->minDelta ))
				{
					listRightBrkPointIntrId++;
				}
				else
				{
					copyElBrkPointIntr(tempListRightBrkPointIntrId, listRightBrkPointIntrId);
					//decreaseByOneRightBrkPointIntr_RR(tempListRightBrkPointIntrId);
					reevaluate_rightBrkPoint_RR(tempListRightBrkPointIntrId, brkPointLeft);
					listRightBrkPointIntrId++;
					tempListRightBrkPointIntrId++;
				}
			}
		}
	}


	ptrMappingOnGenome=g_genomeIndexStart[brkPointLeft];
	while(ptrMappingOnGenome!=NULL)
	{
		newElAdded=true;

		if (ptrMappingOnGenome->readMappingPtr->orientationLeft=='F' && ptrMappingOnGenome->readMappingPtr->orientationRight=='F')
		{
			locBrkPointRightTemp=max(ptrMappingOnGenome->readMappingPtr->libInfo->maxDelta-(brkPointLeft-ptrMappingOnGenome->readMappingPtr->locMapLeftEnd),0)+ptrMappingOnGenome->readMappingPtr->locMapRightEnd;
			locBrkPointLeftTemp =max(ptrMappingOnGenome->readMappingPtr->libInfo->minDelta-(brkPointLeft-ptrMappingOnGenome->readMappingPtr->locMapLeftEnd),0)+ptrMappingOnGenome->readMappingPtr->locMapRightEnd;
		//	printf("L214 %i %i\n",locBrkPointRightTemp, locBrkPointLeftTemp);
		} else if (ptrMappingOnGenome->readMappingPtr->orientationLeft=='R' && ptrMappingOnGenome->readMappingPtr->orientationRight=='R')
		{
			locBrkPointRightTemp=max(ptrMappingOnGenome->readMappingPtr->locMapRightStart-max(ptrMappingOnGenome->readMappingPtr->libInfo->minDelta-(ptrMappingOnGenome->readMappingPtr->locMapLeftStart-brkPointLeft),0),ptrMappingOnGenome->readMappingPtr->locMapLeftStart);
			locBrkPointLeftTemp=max(ptrMappingOnGenome->readMappingPtr->locMapRightStart-max(ptrMappingOnGenome->readMappingPtr->libInfo->maxDelta-(ptrMappingOnGenome->readMappingPtr->locMapLeftStart-brkPointLeft),0),ptrMappingOnGenome->readMappingPtr->locMapLeftStart);
		//	printf("L219 %i %i %i \n", brkPointLeft, locBrkPointRightTemp, locBrkPointLeftTemp);
		}

		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].locBrkPointRight=locBrkPointRightTemp;
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].locBrkPointLeft=locBrkPointLeftTemp;
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].key=locBrkPointLeftTemp;		
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].keyLorR='L';
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].readMappingPtr=ptrMappingOnGenome->readMappingPtr;
		tempListRightBrkPointIntrId++;
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].locBrkPointRight=locBrkPointRightTemp;
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].locBrkPointLeft=locBrkPointLeftTemp;
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].key=locBrkPointRightTemp;
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].keyLorR='R';
		g_tempListRightBrkPointIntr[tempListRightBrkPointIntrId].readMappingPtr=ptrMappingOnGenome->readMappingPtr;
		tempListRightBrkPointIntrId++;
		ptrMappingOnGenome=ptrMappingOnGenome->next;
	}
	temp=g_listRightBrkPointIntr;
	g_listRightBrkPointIntr=g_tempListRightBrkPointIntr;
	g_tempListRightBrkPointIntr=temp;
	g_listRightBrkPointIntrCount=tempListRightBrkPointIntrId;
	return (newElAdded);
}






// This is the head function called for creating clusters. It startes the loop which iterates the left breakpoint
int createInversionClusters(int chroSize) 
{
	int leftBreakPoint; // The value of the left breakpoint considered
	// The genome starts from 0
	int startLeftWindow, endLeftWindow;
	g_listRightBrkPointIntrCount = 0;
	bool newElAdded = false;
	// The loop which iterates on the chromosome repersenting the left breakpoint
	for (leftBreakPoint = 1; leftBreakPoint < chroSize; leftBreakPoint++) 
	{
		newElAdded=false;
		newElAdded = createBreakPointIntervals_Inversion(leftBreakPoint);		
		//if (newElAdded) For Inversions we need to check for every new BrkPoint
		createIntersectingIntervals(leftBreakPoint, 3);
	}
	flushOut(g_listPotClusterFound, leftBreakPoint,3);
}



