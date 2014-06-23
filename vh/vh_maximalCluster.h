#ifndef __VH_CREATE_CLUSTER_DELETION__H__
#define __VH_CREATE_CLUSTER_DELETION__H__
#include "vh_divethandler.h"
#include "vh_gaphandler.h"
#include "vh_heap.h"




typedef struct MappingOnGenome
{
  //int readMappingId;    // the id of the readMapping (starting from 0)
  DivetRow *readMappingPtr;	//
  struct MappingOnGenome *next;
} MappingOnGenome;


typedef struct RightBrkPointInterval
{
  //int locMapLeft;               //the mapping location of the left had
  //int libDeltaMax;      //the loci on the genome which is the starting point of including that mapping in the cluster. locMapLeft-\deltaMax
  int locBrkPointLeft;		// the left end of the right breakpoint interval
  int locBrkPointRight;		// the right end of the right breakpoint interval
  int key;			// Key is the locBrkPointLeft or locBrkPointRight 
  ////////////////////////////////////New field added
  char keyLorR;			// this field indicates that the key is the locBrkPointLeft or locBrkPointRight
  // the values it can get is 'L' or 'R'
  ////////////////////////////////////
  //int readMappingId;
  struct DivetRow *readMappingPtr;
  //      LibraryInfo *libInfo;
} RightBrkPointInterval;


typedef struct ClustersFound
{
  //int readMappingIdArray[MAX_CLUSTER_SIZE]; 
  //DivetRow* readMappingPtrArray[MAX_CLUSTER_SIZE];
  int *readMappingIdArray;
  struct DivetRow **readMappingPtrArray;
  int clusterSize;
  int leftBrkPoint;
  int isMaximalCluster;
  struct ClustersFound *next;
} ClustersFound;


void vh_initializeReadMapping_Deletion (char *, int);
void vh_initializeReadMapping_Inversion (char *, int);
void vh_initializeReadMapping_Insertion (char *, int);
void vh_finalizeReadMapping (char *, int);
int vh_createDeletionClusters (int);
int vh_createInversionClusters (int);
int vh_createInsertionClusters (int);
int vh_compare (const void *, const void *);
int vh_compareReadName (const void *, const void *);
int vh_compareInt (const void *, const void *);
int vh_comparePtr (const void *, const void *);
int vh_noGap (char *, int, int);
int vh_max (int, int);
int vh_min (int, int);
void vh_freeLinkedList (struct MappingOnGenome *);
void vh_finalizeReadMapping (char *, int);
int vh_copyElBrkPointIntr (int, int);
int vh_isItSubset (int *, int, int *, int);
int vh_outputCluster (struct ClustersFound *, int);
int vh_addToPotentialOutput (int, struct Heap *, int);
int vh_flushOut (struct ClustersFound *, int, int);
int vh_createIntersectingIntervals (int, int);
int vh_notBothDirections (struct ClustersFound *);

extern struct MappingOnGenome **g_genomeIndexStart;
extern struct MappingOnGenome **g_genomeIndexEnd;
extern struct RightBrkPointInterval *g_listRightBrkPointIntr;
extern struct RightBrkPointInterval *g_tempListRightBrkPointIntr;
extern int g_maxListBrkPointIntr;
extern int g_listRightBrkPointIntrCount;
extern int g_maxDeltaAmongLibs;
extern struct ClustersFound *g_listPotClusterFound;

#endif
