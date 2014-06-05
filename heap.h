#ifndef _VH_HEAP__H_
#define _VH_HEAP__H_

#define MAX_CLUSTER_SIZE 10000000 //TODO: Ye karish bokon

#include "vh_common.h"
#include <stdbool.h>

typedef struct HeapEl{	
	//int readMappingId;
	DivetRow *readMappingPtr;
	int priorityValue;
} HeapEl;

typedef struct Heap{
	HeapEl heapArray[MAX_CLUSTER_SIZE];
	int heapSize;
} Heap;

extern Heap* g_intersectInterval;

int copyHeapEl(HeapEl *dest, HeapEl *src);
int push_heap(Heap *heapName, HeapEl *newEl);
int minValue_heap(Heap *heapName);
int heapBubleDown(Heap *heapName);
int heap_remove_top(Heap *heapName);
int addToHeap(DivetRow* readMappingPtr, int priorityValue, Heap *heapName);
int writeHeap(Heap *heapName);

#endif

