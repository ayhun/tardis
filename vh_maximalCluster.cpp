#include "vh_maximalCluster.h"
#include "vh_common.h"

//RightBrkPointInterval* g_listRightBrkPointIntr;
MappingOnGenome** g_genomeIndexStart;
MappingOnGenome** g_genomeIndexEnd;
RightBrkPointInterval* g_listRightBrkPointIntr;
RightBrkPointInterval* g_tempListRightBrkPointIntr;
int g_maxListBrkPointIntr;
int g_listRightBrkPointIntrCount;
int g_maxDeltaAmongLibs=0;
ClustersFound *g_listPotClusterFound;
int counter=0;

int min (int x, int y)
{
	if (x<y) return x;
	else return y;
}


int max(int x, int y)
{
	if (x<y) return y;
	else return x;
}

int compare(const void *a, const void *b)
{
	
	struct rightBrkPointInterval *arg1=(struct rightBrkPointInterval *)a;
	struct rightBrkPointInterval *arg2=(struct rightBrkPointInterval *)b;
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

	if ((*arg1).key==(*arg2).key)
	{
		if ((*arg1).keyLorR=='L')
			return -1;
		else return 1;
	}	
}

int compareReadName(const void *a, const void *b)
{
	return strcmp((*(DivetRow **)a)->readName->readName,(*(DivetRow **)b)->readName->readName);
}

int compareInt(const void *a, const void *b)
{
	return *(int*) a - *(int*)b;
}

int comparePtr(const void *a, const void *b)
{
	 
	if (strcmp((*(DivetRow **)a)->readName->readName,(*(DivetRow **)b)->readName->readName)==0)
		{
			if ((*(DivetRow **)a)->locMapLeftEnd==(*(DivetRow **)b)->locMapLeftEnd)
				{
					return ((*(DivetRow **)a)->locMapRightStart - (*(DivetRow **)b)->locMapRightStart);
				}
			else return (*(DivetRow **)a)->locMapLeftEnd - (*(DivetRow **)b)->locMapLeftEnd;
		}
	else return strcmp((*(DivetRow **)a)->readName->readName,(*(DivetRow **)b)->readName->readName); 
}



bool noGap(char *chroName, int startMapping, int endMapping)
{

	for (int count=0; count<g_gapTableSize; count++)
	{
		if (strcmp(chroName, g_gapTable[count].chroName)==0)
		{
			if (startMapping<g_gapTable[count].start && endMapping>g_gapTable[count].end)
				return false;
			if (startMapping>g_gapTable[count].start && startMapping<g_gapTable[count].end)
				return false;
			if (endMapping>g_gapTable[count].start && endMapping<g_gapTable[count].end)
				return false;
		}
	}
	return true;
}


void freeLinkedList(MappingOnGenome* cur)
{
	MappingOnGenome *next;
	while (cur != NULL)
	{
		next = cur->next;
		free(cur);
		cur = next;
	}
}

void finalizeReadMapping(char* chroName, int chroSize)
{
	//Free g_genomeIndexStart and g_genomeIndexEnd and their linked list
	for (int i = 0; i < chroSize; i++)
	{
		freeLinkedList(g_genomeIndexEnd[i]);
		g_genomeIndexEnd[i] = NULL;
		freeLinkedList(g_genomeIndexStart[i]);
		g_genomeIndexStart[i] = NULL;
	}
	free(g_genomeIndexStart);
	g_genomeIndexStart = NULL;
	free(g_genomeIndexEnd);
	g_genomeIndexEnd = NULL;

	//Free g_listRightBrkPointIntrCount and set the g_listRightBrkPointIntrCount
	free(g_listRightBrkPointIntr);
	g_listRightBrkPointIntr = NULL;
	g_listRightBrkPointIntrCount = 0;

	//Free g_tempListRightBrkPointIntr 
	free(g_tempListRightBrkPointIntr);
	g_tempListRightBrkPointIntr = NULL;

	//Free g_intersectInterval -> Heap AND set the heapsize
	free(g_intersectInterval);
	g_intersectInterval = NULL;
}

int copyElBrkPointIntr(int dest, int src)
{
	g_tempListRightBrkPointIntr[dest].locBrkPointLeft=g_listRightBrkPointIntr[src].locBrkPointLeft;
	g_tempListRightBrkPointIntr[dest].locBrkPointRight=g_listRightBrkPointIntr[src].locBrkPointRight;
	g_tempListRightBrkPointIntr[dest].key=g_listRightBrkPointIntr[src].key;
	g_tempListRightBrkPointIntr[dest].readMappingPtr=g_listRightBrkPointIntr[src].readMappingPtr;
	g_tempListRightBrkPointIntr[dest].keyLorR = g_listRightBrkPointIntr[src].keyLorR;
}

bool isItSubset(int *querySet, int querySetSize, int *patternSet, int patternSetSize)
{
	int idPatternSet=0;
	int idQuerySet=0;
	if (patternSetSize==0) return false;
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
			return false;
	}

if (idPatternSet==patternSetSize && idQuerySet<querySetSize)
	return false;
else return true;
}

int flushOut(ClustersFound *listPotClustersFound, int leftBreakPoint, int SVtype)
{
	//printf("%i\n", SVtype);	
	ClustersFound *ptrToOldClusterList;
	ClustersFound *tempPtrOldCluster;
	ptrToOldClusterList=g_listPotClusterFound;
	while (ptrToOldClusterList!=NULL)
	{
		if (ptrToOldClusterList->next!=NULL && ptrToOldClusterList->next->isMaximalCluster==false)
		{
	//		printf("L174\n");
			tempPtrOldCluster=ptrToOldClusterList->next->next;
			free(ptrToOldClusterList->next->readMappingPtrArray);
			free(ptrToOldClusterList->next->readMappingIdArray);
			free(ptrToOldClusterList->next);
			ptrToOldClusterList->next=tempPtrOldCluster;
		}else if (ptrToOldClusterList->next!=NULL && leftBreakPoint > ptrToOldClusterList->next->leftBrkPoint+g_maxDeltaAmongLibs){
	//		printf("L181\n");
			outputCluster(ptrToOldClusterList->next, SVtype);
			tempPtrOldCluster=ptrToOldClusterList->next->next;
			free(ptrToOldClusterList->next->readMappingPtrArray);
			free(ptrToOldClusterList->next->readMappingIdArray);
			free(ptrToOldClusterList->next);
			ptrToOldClusterList->next=tempPtrOldCluster;
		}else
			ptrToOldClusterList=ptrToOldClusterList->next;
	}
	if (g_listPotClusterFound!=NULL && g_listPotClusterFound->isMaximalCluster==false)
	{
		tempPtrOldCluster=g_listPotClusterFound->next;
		free(g_listPotClusterFound->readMappingIdArray);
		free(g_listPotClusterFound->readMappingPtrArray);
		free(g_listPotClusterFound);
		g_listPotClusterFound=tempPtrOldCluster;
	}else if (g_listPotClusterFound!=NULL && leftBreakPoint>g_listPotClusterFound->leftBrkPoint+g_maxDeltaAmongLibs)
		{	
			outputCluster(g_listPotClusterFound, SVtype);
			tempPtrOldCluster=g_listPotClusterFound->next;
			free(g_listPotClusterFound->readMappingPtrArray);
			free(g_listPotClusterFound->readMappingIdArray);
			free(g_listPotClusterFound);
			g_listPotClusterFound=tempPtrOldCluster;
		}

}

int addToPotentialOutput(int leftBreakPoint, Heap *heapName, int SVtype)
{
	
	ClustersFound *newCluster;
	ClustersFound *ptrToOldClusterList;
	ClustersFound *tempPtrOldCluster;
	bool newClusterIsMaximal=true;
	bool oldClusterIsMaximal=true;
	newCluster = (ClustersFound *) malloc(sizeof(ClustersFound));
	newCluster->next=NULL;
	newCluster->isMaximalCluster=true;
	newCluster->leftBrkPoint=leftBreakPoint;
	newCluster->clusterSize=heapName->heapSize;
	newCluster->readMappingPtrArray = (DivetRow **) malloc(newCluster->clusterSize*sizeof(DivetRow*));
	newCluster->readMappingIdArray = (int*) malloc (newCluster->clusterSize * sizeof(int));
	for (int count=0; count<heapName->heapSize; count++)
	{
		newCluster->readMappingPtrArray[count]=heapName->heapArray[count].readMappingPtr;
		newCluster->readMappingIdArray[count]=heapName->heapArray[count].readMappingPtr->divetRowId;
	}
	qsort(newCluster->readMappingIdArray, newCluster->clusterSize, sizeof(int), compareInt);
	ptrToOldClusterList=g_listPotClusterFound;
	while(ptrToOldClusterList!=NULL && newClusterIsMaximal==true)
	{
		if (newClusterIsMaximal==true)
			newClusterIsMaximal=(!isItSubset(newCluster->readMappingIdArray, newCluster->clusterSize, ptrToOldClusterList->readMappingIdArray, ptrToOldClusterList->clusterSize));
		if (newClusterIsMaximal==false)
			ptrToOldClusterList->leftBrkPoint=leftBreakPoint;
		oldClusterIsMaximal=(!isItSubset(ptrToOldClusterList->readMappingIdArray, ptrToOldClusterList->clusterSize, newCluster->readMappingIdArray, newCluster->clusterSize));
		if (oldClusterIsMaximal==false && newClusterIsMaximal==true)
			ptrToOldClusterList->isMaximalCluster=false;
		ptrToOldClusterList=ptrToOldClusterList->next;
	}	


	if (newClusterIsMaximal==true)
	{
		newCluster->next=g_listPotClusterFound;
		g_listPotClusterFound=newCluster;
		ptrToOldClusterList=g_listPotClusterFound;
	}
	if (newClusterIsMaximal==false)
		{
			free(newCluster->readMappingIdArray);
			free(newCluster->readMappingPtrArray);
			free(newCluster);
			newCluster = NULL;
		}

	flushOut(g_listPotClusterFound,leftBreakPoint,SVtype);
}

int outputCluster(ClustersFound *cluster, int SVtype)
{

int listOfReadsOutputed[10000][3];// [0]:locMapLeftEnd, [1]: locMapRightStart, [2]: edit distance. This is used to remove the duplicated reads (clonal) as a result of PCR duplication from each cluster.
int totalAddedToList=0;
bool clonalRead=false;
//	printf("Cluster size %i\n", cluster->clusterSize);
	if (cluster->clusterSize<2)
		return 0;
//		printf("Here\n");
	if (SVtype==3 && notBothDirections(cluster))
		return 0;
	qsort(cluster->readMappingPtrArray, cluster->clusterSize, sizeof(DivetRow **), compareReadName);
//	printf("L229\n");
	for (int readMapCount=0; readMapCount < cluster->clusterSize; readMapCount++)
	{
	clonalRead=false;
		if (readMapCount==0 || strcmp(cluster->readMappingPtrArray[readMapCount]->readName->readName , cluster->readMappingPtrArray[readMapCount-1]->readName->readName)!=0)
			{
			for (int countListOutputed=0; countListOutputed<totalAddedToList; countListOutputed++)
			{
				if (cluster->readMappingPtrArray[readMapCount]->locMapLeftEnd == listOfReadsOutputed[countListOutputed][0] && 
				cluster->readMappingPtrArray[readMapCount]->locMapRightStart == listOfReadsOutputed[countListOutputed][1] &&
				cluster->readMappingPtrArray[readMapCount]->editDistance == listOfReadsOutputed[countListOutputed][2])
					clonalRead=true; 
			}
				if (clonalRead==false)
				{
					fprintf(fileOutput,"%s %s %i %s %i %i %g %g %s %c %c ", cluster->readMappingPtrArray[readMapCount]->readName->readName, cluster->readMappingPtrArray[readMapCount]->chroName, cluster->readMappingPtrArray[readMapCount]->locMapLeftEnd, cluster->readMappingPtrArray[readMapCount]->chroName, cluster->readMappingPtrArray[readMapCount]->locMapRightStart, SVtype, cluster->readMappingPtrArray[readMapCount]->phredScore, cluster->readMappingPtrArray[readMapCount]->editDistance ,cluster->readMappingPtrArray[readMapCount]->libInfo->libName, cluster->readMappingPtrArray[readMapCount]->orientationLeft, cluster->readMappingPtrArray[readMapCount]->orientationRight);
					listOfReadsOutputed[totalAddedToList][0] = cluster->readMappingPtrArray[readMapCount]->locMapLeftEnd;
					listOfReadsOutputed[totalAddedToList][1] = cluster->readMappingPtrArray[readMapCount]->locMapRightStart;
					listOfReadsOutputed[totalAddedToList][2] = (int) cluster->readMappingPtrArray[readMapCount]->editDistance;
					totalAddedToList++;

				}
			}
	}
	//printf("Done\n");
	fprintf(fileOutput,"END\n");
}

int createIntersectingIntervals(int leftBreakPoint, int SVtype)
{
//	printf("l:269\n");
	bool newElAdded;
	newElAdded=false;
	qsort(g_listRightBrkPointIntr, g_listRightBrkPointIntrCount, sizeof(rightBrkPointInterval), compare);
	for (int countIntrEndPoints=0; countIntrEndPoints<g_listRightBrkPointIntrCount; countIntrEndPoints++)
	{
		counter++;
	//	if (counter%10000==0)
//		printf("L302 %i %i %i\n", g_listRightBrkPointIntrCount, g_intersectInterval->heapSize, leftBreakPoint);
		
		//if (g_listRightBrkPointIntr[countIntrEndPoints].locBrkPointLeft==g_listRightBrkPointIntr[countIntrEndPoints].key)
		if (g_listRightBrkPointIntr[countIntrEndPoints].keyLorR=='L')
		{
			newElAdded=true;
			addToHeap(g_listRightBrkPointIntr[countIntrEndPoints].readMappingPtr, g_listRightBrkPointIntr[countIntrEndPoints].locBrkPointRight, g_intersectInterval);
			
//if (leftBreakPoint>500 && leftBreakPoint<1500)
//	printf("%i %i\n", leftBreakPoint, g_intersectInterval->heapSize);

		//for (int testCount=0; testCount<g_intersectInterval->heapSize; testCount++)
		//{
		//}
		//writeHeap(g_intersectInterval);

		} else if (g_listRightBrkPointIntr[countIntrEndPoints].keyLorR=='R')
			//(g_listRightBrkPointIntr[countIntrEndPoints].locBrkPointRight==g_listRightBrkPointIntr[countIntrEndPoints].key)
		{
			if (minValue_heap(g_intersectInterval)==g_listRightBrkPointIntr[countIntrEndPoints].key)
			{
				if (newElAdded==true)
				{
					addToPotentialOutput(leftBreakPoint, g_intersectInterval, SVtype);	
					newElAdded=false;
				}
				//writeHeap(g_intersectInterval);
				heap_remove_top(g_intersectInterval);
			}else 
				{
					printf("An Error Occured\n");
					writeHeap(g_intersectInterval);
				}
		}
	}
/*if (leftBreakPoint>500 && leftBreakPoint<1500)
	printf("%i %i\n", leftBreakPoint, g_intersectInterval->heapSize);
*/
	g_intersectInterval->heapSize=0;
}

bool notBothDirections(ClustersFound *cluster)
{
	bool FF,RR;
	FF=false;
	RR=false;
	for(int readMapCount=0; readMapCount < cluster->clusterSize; readMapCount++)
	{
		if (cluster->readMappingPtrArray[readMapCount]->orientationLeft=='F' && cluster->readMappingPtrArray[readMapCount]->orientationRight=='F')
			FF=true;

		if (cluster->readMappingPtrArray[readMapCount]->orientationLeft=='R' && cluster->readMappingPtrArray[readMapCount]->orientationRight=='R')
			RR=true;

	if (FF && RR)
		return false;
	}
return true;
}

