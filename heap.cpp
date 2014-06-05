#include "vh_heap.h"
#include "vh_divethandler.h"
Heap *g_intersectInterval;



int writeHeap(Heap *heapName)
{
	for(int count=0; count<heapName->heapSize; count++)
	{
		printf("%i %i %i *", heapName->heapArray[count].priorityValue, heapName->heapArray[count].readMappingPtr->locMapLeftEnd, heapName->heapArray[count].readMappingPtr->locMapRightStart);
	}
	printf("END\n");
}


int copyHeapEl(HeapEl *dest, HeapEl *src)
{
	dest->readMappingPtr = src->readMappingPtr;
	dest->priorityValue = src->priorityValue;
}


int push_heap(Heap *heapName, HeapEl *newEl)
{
	int heapIndex = heapName->heapSize;
	HeapEl tempEl;
	heapName->heapSize++;
	copyHeapEl(&(heapName->heapArray[heapIndex]), newEl);
	heapIndex = heapName->heapSize - 1;;
	while(heapIndex > 0 && heapName->heapArray[heapIndex].priorityValue < heapName->heapArray[(heapIndex + 1) / 2 -1].priorityValue)
	{
		copyHeapEl(&tempEl, &heapName->heapArray[heapIndex]);
		copyHeapEl(&(heapName->heapArray[heapIndex]),  &(heapName->heapArray[(heapIndex + 1) / 2 -1]));
		copyHeapEl(&(heapName->heapArray[(heapIndex + 1) / 2 -1]), &tempEl);
		heapIndex = (heapIndex + 1) / 2 - 1;
	}
}

int minValue_heap(Heap *heapName)
{
	if (heapName->heapSize > 0)
		return (heapName->heapArray[0].priorityValue);
	else 
		return -1000;
}

int heapBubleDown(Heap *heapName)
{
	int heapIndex = 0;
	HeapEl tempEl;
	while((heapIndex * 2 + 1<heapName->heapSize && heapName->heapArray[heapIndex].priorityValue > heapName->heapArray[heapIndex * 2 + 1].priorityValue)||(heapIndex * 2 + 2 < heapName->heapSize && heapName->heapArray[heapIndex].priorityValue > heapName->heapArray[heapIndex * 2 + 2].priorityValue))
	{

		if (heapIndex * 2 + 2<heapName->heapSize)
		{
			if (heapName->heapArray[heapIndex * 2 + 1].priorityValue > heapName->heapArray[heapIndex * 2 + 2].priorityValue)
			{
				copyHeapEl(&tempEl, & (heapName->heapArray[heapIndex]));
				copyHeapEl(& (heapName->heapArray[heapIndex]), &(heapName->heapArray[heapIndex * 2 + 2]));
				copyHeapEl(& (heapName->heapArray[heapIndex * 2 + 2]), &tempEl);
				heapIndex = heapIndex * 2 + 2;

			}
			else
			{
				copyHeapEl(&tempEl, & (heapName->heapArray[heapIndex]));
				copyHeapEl(& (heapName->heapArray[heapIndex]), &(heapName->heapArray[heapIndex * 2 + 1]));
				copyHeapEl(&(heapName->heapArray[heapIndex * 2 + 1]), &tempEl);
				heapIndex = heapIndex * 2 + 1;
			}
		} 
		else
		{
			copyHeapEl(&tempEl, & (heapName->heapArray[heapIndex]));
			copyHeapEl(& (heapName->heapArray[heapIndex]), &(heapName->heapArray[heapIndex * 2 + 1]));
			copyHeapEl(&(heapName->heapArray[heapIndex * 2 + 1]), &tempEl);
			heapIndex = heapIndex * 2 + 1;
		}
	}
}

int heap_remove_top(Heap *heapName)
{


	if (heapName->heapSize > 0)
	{
		copyHeapEl(&(heapName->heapArray[0]), &(heapName->heapArray[heapName->heapSize - 1]));
		heapBubleDown(heapName);
		heapName->heapSize--;
	}
}

int addToHeap(DivetRow* readMappingPtr, int priorityValue, Heap *heapName)
{
	HeapEl *newEl;
	newEl = (HeapEl *)malloc(sizeof(HeapEl));
	newEl->readMappingPtr = readMappingPtr;
	newEl->priorityValue = priorityValue;
	push_heap(heapName, newEl);
	free(newEl);
}

