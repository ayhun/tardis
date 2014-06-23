#include <stdio.h>
#include <stdlib.h>
#include <string.h>


const float weight_1_2_3 = 5;
const float weight_1_2 = 61;
const float weight_1_3 = 10;
const float weight_2_3 = 9;
const float weight_1 = 14;
const float weight_2 = 12;
const float weight_3 = 400;

const int maxLengthSV = 500000;
const int maxListClusterSize = 3000000;
const int strSize = 100;
const int maxNumSV = 1000000;

char multiInd[10][strSize];
int multiIndCount;

typedef struct multiLib
{
  char libName[strSize];
  int indId;			// index to the multiInd[X] 
  int maxInstSize;
  int minInstSize;
  int readLen;
} multiLib;

multiLib multiLibs[10];
int multiLibsCount;

typedef struct SV_selected
{
  char chroName[strSize];
  char SVtype;			// D: deletion, V: Inversion, I: insertion
  int posStart_SV, posEnd_SV;
  int sup[10];
  SV_selected *conflict_Next;
} SV_selected;


SV_selected listSelectedSV[maxNumSV];
int numSV;

////////////////////////////////////////////////


typedef struct clusterIdEl
{
  int clusterId;
  clusterIdEl *next;
};




typedef struct readEl
{
  char readName[strSize];
  int readId;
  int indId;
  int readCovered;		//0 not covered, 1 covered;
  clusterIdEl *next;
};



readEl *listReadEl;
int sizeListReadEl;

///////////////////////////////////////////////


typedef struct readMappingEl
{
  int readId;
  char chroName[strSize];
  int posMapLeft;
  int posMapRight;
  char orient1;
  char orient2;
  int indId;
  float probEditDist;
  readMappingEl *next;
};




typedef struct clusterEl
{
  int clusterId;
  char chroName[strSize];
  int posStartSV;
  int posEndSV;
  char SVtype;			//V: inversion, D: Deletion, I: insertion
  int picked;			// This nodes has been picked
  //int canNotBePicked;// Two nodes in a conflict triangle have been picked
  //int indIdCount[10];// count number of reads from each individuals you have
  int indIdCount[10];
  readMappingEl *next;
//      clusterIdEl *nextGraphNode; // the pointer for the conflict graph

};

clusterEl *listClusterEl;
int sizeListClusterEl;

//////////////////////////////////////////
int numCallsRequsted;

int
max (int a, int b)
{
  if (a < b)
    return b;
  else
    return a;
}


int
min (int a, int b)
{
  if (a < b)
    return a;
  else
    return b;
}

int vh_processTheSV (int clusterId)
{
  readMappingEl *readMappingElPtr;
  readMappingElPtr = listClusterEl[clusterId].next;
  int posStartSV = 0, posEndSV = 0;
  int posStartSV2 = 300000000, posEndSV2 = 300000000;
  int maxDelLen, minDelLen;
  strcpy (listClusterEl[clusterId].chroName, readMappingElPtr->chroName);
  maxDelLen = 1000000000;
  minDelLen = 0;
  if (listClusterEl[clusterId].SVtype == 'V')
    {
      while (readMappingElPtr != NULL)
	{

	  if (readMappingElPtr->posMapLeft > posStartSV
	      && readMappingElPtr->orient1 == 'F')
	    posStartSV = readMappingElPtr->posMapLeft;

	  if (readMappingElPtr->posMapLeft < posStartSV2
	      && readMappingElPtr->orient1 == 'R')
	    posStartSV2 = readMappingElPtr->posMapLeft;

	  if (readMappingElPtr->posMapRight > posEndSV
	      && readMappingElPtr->orient1 == 'F')
	    posEndSV = readMappingElPtr->posMapRight;

	  if (readMappingElPtr->posMapRight < posEndSV2
	      && readMappingElPtr->orient1 == 'R')
	    posEndSV2 = readMappingElPtr->posMapRight;
	  readMappingElPtr = readMappingElPtr->next;

	}
      listClusterEl[clusterId].posStartSV = posStartSV;
      listClusterEl[clusterId].posEndSV = posEndSV;
      for (int brk1 = posStartSV; brk1 <= posStartSV2; brk1++)
	{
	  for (int brk2 = posEndSV; brk2 <= posEndSV2; brk2++)
	    {
	      readMappingElPtr = listClusterEl[clusterId].next;
	      while (readMappingElPtr != NULL)
		{
		  if ((readMappingElPtr->orient1 == 'F'
		       && brk2 - readMappingElPtr->posMapRight + brk1 -
		       readMappingElPtr->posMapLeft > 20
		       && brk2 - readMappingElPtr->posMapRight + brk1 -
		       readMappingElPtr->posMapLeft < 113)
		      || (readMappingElPtr->orient1 == 'R'
			  && 0 - brk2 + readMappingElPtr->posMapRight - brk1 +
			  readMappingElPtr->posMapLeft > 20
			  && 0 - brk2 + readMappingElPtr->posMapRight - brk1 +
			  readMappingElPtr->posMapLeft < 113))
		    {

		      readMappingElPtr = readMappingElPtr->next;
		      if (readMappingElPtr == NULL)
			{
			  //              printf("%i %i %i %i %i %i\n", brk1, brk2, posStartSV, posStartSV2, posEndSV, posEndSV2);
			  return 0;
			}

		    }
		  else
		    readMappingElPtr = NULL;

		}
	    }
	}

      printf ("ERRORORORORORO %i %i %i %i %i\n", posStartSV, posStartSV2,
	      posEndSV, posEndSV2, clusterId);
    }
}






int vh_init (FILE * fpCluster)
{
  char readNameStr[strSize];
  char chroName[strSize];
  int startPos;
  int stopPos;
  int SVtype;
  float editProbDist;
  char libName[strSize], indName[strSize], filePath[strSize];
  int listReadElId = 0, listClusterElId = 0;
  int readId, minInstSize, maxInstSize, readLen;
  char orient1, orient2;
  clusterIdEl *clusterIdElNew;
  readMappingEl *readMappingElNew;
  multiLibsCount = 0;
  int numLibInFile;
  listClusterEl =
    (clusterEl *) malloc (maxListClusterSize * sizeof (clusterEl));
  for (int count = 0; count < maxListClusterSize; count++)
    {
      listClusterEl[count].clusterId = 0;
      listClusterEl[count].next = NULL;
      //listClusterEl[count].nextGraphNode=NULL;
    }
  printf ("L264\n");
  listClusterElId = 0;

  listClusterEl[listClusterElId].clusterId = listClusterElId;
  listClusterEl[listClusterElId].picked = 0;
  //listClusterEl[listClusterElId].canNotBePicked=0;



  while (fscanf (fpCluster, "%s ", readNameStr) != EOF)
    {
      if (strcmp (readNameStr, "END") == 0)
	{
	  //      processTheSV(listClusterElId);
	  //      addToConflictGraph(listClusterElId);
	  /*listClusterElId++;
	     if ((listClusterElId%10000)==0)
	     printf("%i\n", listClusterElId);
	     listClusterEl[listClusterElId].clusterId=listClusterElId;
	     listClusterEl[listClusterElId].picked=0; */
	  if (SVtype == 3)
	    listClusterEl[listClusterElId].SVtype = 'V';
	  if (SVtype == 2)
	    listClusterEl[listClusterElId].SVtype = 'D';
	  if (SVtype == 1)
	    listClusterEl[listClusterElId].SVtype = 'I';
	  vh_processTheSV (listClusterElId);
	  listClusterElId++;
	  if ((listClusterElId % 10000) == 0)
	    printf ("%i\n", listClusterElId);
	  listClusterEl[listClusterElId].clusterId = listClusterElId;
	  listClusterEl[listClusterElId].picked = 0;

	}
      else
	{
	  fscanf (fpCluster, "%s %i %i %i %f %s %c %c ", chroName, &startPos,
		  &stopPos, &SVtype, &editProbDist, libName, &orient1,
		  &orient2);
	  /////////////////////////////ADD THE CLUSTER TO THE READ////////////////////////
	  //readId=findReadName(readNameStr);
	  //listReadEl[readId].indId=findIndId(libName);
	  //clusterIdElNew = (clusterIdEl *) malloc(sizeof(clusterIdEl));
	  //printf("%s\n",        readNameStr);
	  //clusterIdElNew->clusterId=listClusterElId;
	  //clusterIdElNew->next=listReadEl[readId].next;
	  //listReadEl[readId].next=clusterIdElNew;
	  //listReadEl[readId].indId=findIndId(libName);
	  //////////////////////////////ADD THE READ TO THE CLUSTER/////////////////////////
	  readMappingElNew =
	    (readMappingEl *) malloc (sizeof (readMappingEl));
	  //listClusterEl[listClusterElId].clusterId=listClusterElId;
	  //listClusterEl[listClusterElId].picked=0;
	  //listClusterEl[listClusterElId].canNotBePicked=0;
	  //listClusterEl[listClusterElId].indId[findIndId(libName)]=0;
	  //listClusterEl[listClusterElId].indIdCount[findIndId(libName)]++;
	  readMappingElNew->readId = readId;
	  //readMappingElNew->indId=findIndId(libName);
	  readMappingElNew->orient1 = orient1;
	  readMappingElNew->orient2 = orient2;
	  readMappingElNew->probEditDist = editProbDist;
	  readMappingElNew->posMapLeft = startPos;
	  readMappingElNew->posMapRight = stopPos;
	  strcpy (readMappingElNew->chroName, chroName);
	  readMappingElNew->next = listClusterEl[listClusterElId].next;
	  listClusterEl[listClusterElId].next = readMappingElNew;
	}
    }

  sizeListClusterEl = listClusterElId;

}



int main (int argv, char *argc[])
{
  FILE *readFp, *clusterFp, *libFp;
  clusterFp = fopen (argc[1], "r");
  vh_init (clusterFp);
//      createConflictGraph();  
//      pickSet();      

/*	for (int count=0; count<sizeListClusterEl; count++)
	{
		
	}
*/


}
