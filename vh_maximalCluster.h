#ifndef __VH_CREATE_CLUSTER_DELETION__H__
#define __VH_CREATE_CLUSTER_DELETION__H__
#include "vh_divethandler.h"
#include "vh_gaphandler.h"
#include "vh_heap.h"
#include <stdbool.h>


const int maxInversionLen=2000000;
const int maxDeletionLen=500000;

typedef struct MappingOnGenome{
	//int readMappingId;	// the id of the readMapping (starting from 0)
	DivetRow *readMappingPtr; //
	MappingOnGenome *next;
} MappingOnGenome;


typedef struct rightBrkPointInterval{
	//int locMapLeft;		//the mapping location of the left had
	//int libDeltaMax;	//the loci on the genome which is the starting point of including that mapping in the cluster. locMapLeft-\deltaMax
	int locBrkPointLeft;	// the left end of the right breakpoint interval
	int locBrkPointRight;	// the right end of the right breakpoint interval
	int key; 		// Key is the locBrkPointLeft or locBrkPointRight 
	////////////////////////////////////New field added
	char keyLorR; // this field indicates that the key is the locBrkPointLeft or locBrkPointRight
				  // the values it can get is 'L' or 'R'
	////////////////////////////////////
	//int readMappingId;
	DivetRow *readMappingPtr;
	//	LibraryInfo *libInfo;
} RightBrkPointInterval;


typedef struct ClustersFound{
	//int readMappingIdArray[MAX_CLUSTER_SIZE]; 
	//DivetRow* readMappingPtrArray[MAX_CLUSTER_SIZE];
	int *readMappingIdArray;
	DivetRow** readMappingPtrArray;
	int clusterSize;
	int leftBrkPoint;
	bool isMaximalCluster;
	ClustersFound *next;
	
}ClustersFound;


void initializeReadMapping_Deletion(char *, int);
void initializeReadMapping_Inversion(char *, int);
void initializeReadMapping_Insertion(char *, int);
void finalizeReadMapping(char *, int);
int createDeletionClusters(int );
int createInversionClusters(int );
int createInsertionClusters(int );
int compare(const void *, const void *);
int compareReadName(const void *, const void *);
int compareInt(const void *, const void *);
int comparePtr(const void *, const void *);
bool noGap(char *, int, int);
int max(int, int);
int min(int, int);
void freeLinkedList(MappingOnGenome *);
void finalizeReadMapping(char *, int);
int copyElBrkPointIntr(int, int);
bool isItSubset(int *, int, int *, int);
int outputCluster(ClustersFound *, int);
int addToPotentialOutput(int, Heap *, int);
int flushOut(ClustersFound *, int, int);
int createIntersectingIntervals(int, int);
bool notBothDirections(ClustersFound*);

extern MappingOnGenome** g_genomeIndexStart;
extern MappingOnGenome** g_genomeIndexEnd;
extern RightBrkPointInterval* g_listRightBrkPointIntr;
extern RightBrkPointInterval* g_tempListRightBrkPointIntr;
extern int g_maxListBrkPointIntr;
extern int g_listRightBrkPointIntrCount;
extern int g_maxDeltaAmongLibs;
extern ClustersFound *g_listPotClusterFound;

#endif

