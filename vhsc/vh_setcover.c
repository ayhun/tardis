#include "vh_setcover.h"
#include "vh_buffer.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

float *weightsForCombination;

int minimumSupNeeded=2;
float mismatchPenalty=0;

char multiInd[totalNumInd][strSize]; //Name of each individual
int multiIndCount;//total number of numtiple individual <= totalNumInd

int coverageForIndividual[totalNumInd];
int totalCoverage;

clusterElRead clusterElRead_Single; // the name of cluster which each new cluster is read to it
clusterEl *listClusterEl; // the array of all the cluster reads
int sizeListClusterEl; // total number of clusters
int conflictResFlag = 0;
int weightedHeuristicFlag = 0;

typedef struct mobileEl{
  int startPos;
  int endPos;
  char nameMobile[50];
}mobileEl;

mobileEl listMobileEl[1000];
int listMobileElSize;



//////////////////////////////////////////
int numCallsRequested; //maximum number of SVs the users wants to us to output

int tempClusterId;

int compareReadMappingEl(const void *a, const void *b) // the compare function for qsort procudre, should sort in increasing order
{
  struct readMappingEl *arg1=(struct readMappingEl *)a;
  struct readMappingEl *arg2=(struct readMappingEl *)b;
  return ((*arg1).editDistance - (*arg2).editDistance);
}

int findLibId(int multiLibsCount, char *libName)
{	
  int count=0;
  for ( count=0; count<multiLibsCount; count++)
    {
      if (strcmp(libName, multiLibs[count].libName)==0)
	return count;
    }
}

readMappingEl *addToListOfReads(readMappingEl *linkList, readMappingEl element) //Adds a nw readMapping to the list of read mapping *linkList
{
  readMappingEl *newEl;
  newEl = (readMappingEl *) malloc(sizeof(readMappingEl));
  newEl->readId = element.readId;
  strcpy(newEl->chroName, element.chroName);
  newEl->posMapLeft = element.posMapLeft;
  newEl->posMapRight = element.posMapRight;
  newEl->indId = element.indId;
  newEl->orient1 = element.orient1;
  newEl->orient2 = element.orient2;
  newEl->editDistance = element.editDistance;
  newEl->probEditDist = element.probEditDist;
  newEl->next=linkList;
  return newEl;
} 



int processTheSV(int clusterId) // each new cluster read (clusterElRead_Single) is copied into a new cell in arra of clusters (listClusterEl) in index (clusterId)
{
  readMappingEl *readMappingElPtr;
  int posStartSV=-1, posEndSV=maxChroSize;
  int posStartSV_Outer=maxChroSize, posEndSV_Outer=-1;
  qsort(clusterElRead_Single.readMappingElArray, clusterElRead_Single.sizeOfCluster, sizeof(readMappingEl), compareReadMappingEl);
  listClusterEl[clusterId].clusterId=clusterId;
  listClusterEl[clusterId].readMappingSelected=NULL;
  strcpy(listClusterEl[clusterId].chroName, clusterElRead_Single.readMappingElArray[0].chroName);
  listClusterEl[clusterId].next=NULL;
  int count=0;
  for (count=0; count<totalNumInd; count++)
    {
      listClusterEl[clusterId].indIdCount[count]=0;
    }

  if (listClusterEl[clusterId].SVtype=='A')
    {
      int count=0;
		
      for (count=0; count<clusterElRead_Single.sizeOfCluster; count++)
	{
	  //Move each paired-end read in the cluster (clusterElRead_Single) into 
	  //the arry of clusters listClusterEl 
	  if (clusterElRead_Single.readMappingElArray[count].orient1=='F')
	    {
	      if (clusterElRead_Single.readMappingElArray[count].posMapLeft > posStartSV)
		posStartSV = clusterElRead_Single.readMappingElArray[count].posMapLeft;
	      if (clusterElRead_Single.readMappingElArray[count].posMapLeft < posStartSV_Outer)
		posStartSV_Outer = clusterElRead_Single.readMappingElArray[count].posMapLeft;
				
	    }
	  if (clusterElRead_Single.readMappingElArray[count].orient1=='R')
	    {
	      if (clusterElRead_Single.readMappingElArray[count].posMapLeft < posEndSV)
		posEndSV = clusterElRead_Single.readMappingElArray[count].posMapLeft;
	      if (clusterElRead_Single.readMappingElArray[count].posMapLeft > posEndSV_Outer)
		posEndSV_Outer = clusterElRead_Single.readMappingElArray[count].posMapLeft;
	    }
	  listClusterEl[clusterId].next = addToListOfReads(listClusterEl[clusterId].next, clusterElRead_Single.readMappingElArray[count]); 
	  if (listClusterEl[clusterId].next==NULL)
	    printf("L73 Null\n");
	}

      listClusterEl[clusterId].posStartSV = posStartSV;
      listClusterEl[clusterId].posStartSV_Outer = posStartSV_Outer;
      listClusterEl[clusterId].posEndSV = posEndSV;
      listClusterEl[clusterId].posEndSV_Outer = posEndSV_Outer;
    } else if(listClusterEl[clusterId].SVtype=='B')
    {
      int count=0;

      for ( count=0; count<clusterElRead_Single.sizeOfCluster; count++)
	{
	  //Move each paired-end read in the cluster (clusterElRead_Single) into 
	  //the arry of clusters listClusterEl 
	  if (clusterElRead_Single.readMappingElArray[count].orient1=='F')
	    {
	      if (clusterElRead_Single.readMappingElArray[count].posMapLeft > posStartSV)
		posStartSV = clusterElRead_Single.readMappingElArray[count].posMapLeft;
	      if (clusterElRead_Single.readMappingElArray[count].posMapLeft < posStartSV_Outer)
		posStartSV_Outer = clusterElRead_Single.readMappingElArray[count].posMapLeft;
	    }
	  if (clusterElRead_Single.readMappingElArray[count].orient1=='R')
	    {
	      if (clusterElRead_Single.readMappingElArray[count].posMapLeft < posEndSV)
		posEndSV = clusterElRead_Single.readMappingElArray[count].posMapLeft;
	      if (clusterElRead_Single.readMappingElArray[count].posMapLeft > posEndSV_Outer)
		posEndSV_Outer = clusterElRead_Single.readMappingElArray[count].posMapLeft;
	    }
	  listClusterEl[clusterId].next = addToListOfReads(listClusterEl[clusterId].next, clusterElRead_Single.readMappingElArray[count]); 
	  if (listClusterEl[clusterId].next==NULL)
	    printf("L73 Null\n");
	}

      listClusterEl[clusterId].posStartSV = posStartSV;
      listClusterEl[clusterId].posStartSV_Outer = posStartSV_Outer;
      listClusterEl[clusterId].posEndSV = posEndSV;
      listClusterEl[clusterId].posEndSV_Outer = posEndSV_Outer;
    }
  else
    {
      int count=0;
      
      for ( count=0; count<clusterElRead_Single.sizeOfCluster; count++)
	{
	  //Move each paired-end read in the cluster (clusterElRead_Single) into 
	  //the arry of clusters listClusterEl 
	  if (clusterElRead_Single.readMappingElArray[count].posMapLeft > posStartSV)
	    posStartSV = clusterElRead_Single.readMappingElArray[count].posMapLeft;
	  if (clusterElRead_Single.readMappingElArray[count].posMapLeft < posStartSV_Outer)
	    posStartSV_Outer = clusterElRead_Single.readMappingElArray[count].posMapLeft;
	  if (clusterElRead_Single.readMappingElArray[count].posMapRight < posEndSV)
	    posEndSV = clusterElRead_Single.readMappingElArray[count].posMapRight;
	  if (clusterElRead_Single.readMappingElArray[count].posMapRight > posEndSV_Outer)
	    posEndSV_Outer = clusterElRead_Single.readMappingElArray[count].posMapRight;
	  listClusterEl[clusterId].next = addToListOfReads(listClusterEl[clusterId].next, clusterElRead_Single.readMappingElArray[count]); 
	}
      
      listClusterEl[clusterId].posStartSV = posStartSV;
      listClusterEl[clusterId].posStartSV_Outer = posStartSV_Outer;
      listClusterEl[clusterId].posEndSV = posEndSV;
      listClusterEl[clusterId].posEndSV_Outer = posEndSV_Outer;
    }
}

int findReadName(char *readName)
{
  int minId=0;
  int maxId=sizeListReadEl;
  int middleId=sizeListReadEl/2;
  if (strcmp(readName, listReadEl[minId].readName)==0)
    return minId;
  if (strcmp(readName, listReadEl[maxId].readName)==0)
    return maxId;
	
  while (strcmp(readName, listReadEl[middleId].readName)!=0)
    {
      if (strcmp(readName, listReadEl[middleId].readName)>0)
	{
	  minId=middleId;
	  middleId=(minId+maxId)/2;
	}
      else if (strcmp(readName, listReadEl[middleId].readName)<0)
	{
	  maxId=middleId;
	  middleId=(minId+maxId)/2;
	}
      if (strcmp(readName, listReadEl[minId].readName)==0)
	return minId;
      if (strcmp(readName, listReadEl[maxId].readName)==0)
	return maxId;
    }
  return middleId;
}

int findIndId(int multiLibsCount, char *libName)
{
  int count=0;

  for (count=0; count<multiLibsCount; count++)
    {
      if (strcmp(multiLibs[count].libName, libName)==0)
	return (multiLibs[count].indId);
    }
}

int addNewInd(char *indName)
{
  int count=0;

  for (count=0; count<multiIndCount; count++)
    {
      if (strcmp(multiInd[count], indName)==0)
	return count;
    }
  strcpy(multiInd[multiIndCount],indName);
  multiIndCount++;
  return multiIndCount-1;
}


int freeLinkList(readMappingEl *ptr)
{
  if (ptr!=NULL) return 0;
  free(ptr);
}



int init(bam_info **in_bams, int num_bams, FILE *fpRead ,FILE *fpCluster, FILE *fpWeights, FILE *fpCoverage)
{

  char orient1, orient2, readNameStr[strSize], chroName1[strSize], libName[strSize], indName[strSize], filePath[strSize], chroName2[strSize];
  int startPos, stopPos, SVtype, listReadElId=0, listClusterElId=0, readId, minInstSize, maxInstSize, readLen, combinationId, indId, editDist, indCov;
  int i,j; int multiLibsCount;
  float editProbDist, weight;
  clusterIdEl *clusterIdElNew;
  readMappingEl *readMappingElNew;
  readMappingEl *ptrMapping;
  fscanf(fpRead, "%i\n",	&sizeListReadEl);
  listReadEl = (readEl *) malloc(sizeListReadEl*sizeof(readEl));
  
  multiLibsCount=0;

  for (i=0; i< num_bams; i++)
    multiLibsCount = in_bams[i]->num_libraries;

  multiLibs = (multiLib *) malloc(sizeof (multiLib) * multiLibsCount);

  multiLibsCount=0;

  for (j=0; j<num_bams; j++)
    {
      for (i=0; i<in_bams[j]->num_libraries; i++)
	{
	  strcpy(multiLibs[multiLibsCount].libName, in_bams[j]->libraries[i]->libname);
	  multiLibs[multiLibsCount].indId=addNewInd(in_bams[j]->sample_name);    
	  multiLibs[multiLibsCount].maxInstSize = in_bams[j]->libraries[i]->conc_max - 2 * in_bams[j]->libraries[i]->read_length; 
	  multiLibs[multiLibsCount].minInstSize = in_bams[j]->libraries[i]->conc_min - 2 * in_bams[j]->libraries[i]->read_length; 
	  if (multiLibs[multiLibsCount].maxInstSize < 0)
	    multiLibs[multiLibsCount].maxInstSize = 0;
	  if (multiLibs[multiLibsCount].minInstSize < 0)
	    multiLibs[multiLibsCount].minInstSize = 0;
	  multiLibs[multiLibsCount].readLen = in_bams[j]->libraries[i]->read_length; 
	  multiLibsCount++;
	}
    }
  

  weightsForCombination = (float *) malloc( (int) pow(2,multiIndCount)*sizeof(float));
  //	printf("%i \n", multiIndCount);
  int count=0;

  for ( count=0; count<pow(2, multiIndCount); count++)
    {
      weightsForCombination[count]=0;
    }	
  if (fpWeights!=NULL)
    {
      int count=0;

      for(count=0; count<pow(2, multiIndCount)-1; count++)
	{
	  //printf("%i\n", count);
	  combinationId=0;
	  fscanf(fpWeights,"%s ", indName);
	  while(strcmp(indName,":")!=0)
	    {
	      indId=addNewInd(indName);		
	      combinationId = combinationId + (int)pow(2,indId);
	      fscanf(fpWeights,"%s ", indName);
			
	    }
	  fscanf(fpWeights, "%f\n", &weight);
	  weightsForCombination[combinationId]=weight;
	}
    } else {	
    int count=0;

    for ( count=1; count<pow(2, multiIndCount); count++)
      {
	weightsForCombination[count]=1; //THE Defualt weights for SVs for combination of individuals 
      }
  }

  if (fpCoverage!=NULL)
    {
      while(fscanf(fpCoverage,"%s %i\n", indName, &indCov)!=EOF)
	{
	  //		printf("%s %i\n", indName, indCov);
	  int count=0;
	  for (count=0; count<multiIndCount; count++)
	    {
	      if(strcmp(multiInd[count], indName)==0)
		{
		  coverageForIndividual[count] = indCov;
		  totalCoverage = totalCoverage + indCov;
		  //	printf("%i\n", totalCoverage);
		}else ; //printf("%s %s\n", multiInd[count], indName);
	    }		
	}

    }else{
    int count=0;

    for (count=0; count<multiIndCount; count++)
      {
	coverageForIndividual[count]=1;
	//printf("coverage %i\n",coverageForIndividual[count]);
	totalCoverage++;
      }
  }

  weightsForCombination[0]=0;
	
  while(fscanf(fpRead, "%s\n", readNameStr)!=EOF)
    {
      strcpy(listReadEl[listReadElId].readName, readNameStr);
      listReadEl[listReadElId].readCovered=0;
      listReadEl[listReadElId].readId=listReadElId;
      listReadEl[listReadElId].next=NULL;
      listReadElId++;
    }
  listClusterEl = (clusterEl *) malloc(maxClustersAllowed*sizeof(clusterEl));

  for (count=0; count<maxClustersAllowed; count++)
    {
      listClusterEl[count].clusterId=0;
      listClusterEl[count].next=NULL;
    }
  listClusterElId=0;

  listClusterEl[listClusterElId].clusterId=listClusterElId;
  for ( count=0; count < multiIndCount; count++)
    listClusterEl[listClusterElId].indIdCount[count]=0;

	

  while(fscanf(fpCluster, "%s ", readNameStr)!=EOF)
    {
      if (strcmp(readNameStr,"END")==0)
	{
	  if (SVtype==3) 
	    listClusterEl[listClusterElId].SVtype='V';
	  if (SVtype==2)
	    listClusterEl[listClusterElId].SVtype='D';
	  if (SVtype==1)
	    listClusterEl[listClusterElId].SVtype='I';
	  if (SVtype==10) // Insertion from chrN in Forward orientation
	    listClusterEl[listClusterElId].SVtype='A';
	  if (SVtype==12) // Insertion from chrN in Reverse orientation
	    listClusterEl[listClusterElId].SVtype='B';
	  //listClusterEl[listClusterElId].next = (readMappingEl *) malloc(multiIndCount*sizeof(readMappingEl));
	  ///////////////////If a new cluster is read indicates an SV bigger than what is expection (user defined threshold) it is pruned out
			
	  processTheSV(listClusterElId);
	  tempClusterId++;
	  if ((listClusterEl[listClusterElId].posEndSV-listClusterEl[listClusterElId].posStartSV > maxLengthSV_Del && listClusterEl[listClusterElId].SVtype=='D' )
	      ||( listClusterEl[listClusterElId].posEndSV-listClusterEl[listClusterElId].posStartSV > maxLengthSV_Inv && listClusterEl[listClusterElId].SVtype=='V' ))
	    {
	      //TROW THE LINK LIST OUT

	      ptrMapping=listClusterEl[listClusterElId].next;
				
	      while(ptrMapping!=NULL)
		{
		  listReadEl[ptrMapping->readId].next=listReadEl[ptrMapping->readId].next->next;
		  //free(tmpPtrMapping);
		  ptrMapping=ptrMapping->next;
		}
	      ptrMapping=listClusterEl[listClusterElId].next;
	      freeLinkList(ptrMapping);
				
	    } else 
	    {
	      listClusterElId++;
	    }
			
	  listClusterEl[listClusterElId].clusterId=listClusterElId;
	  int multiIndId;
	  for (multiIndId=0; multiIndId<totalNumInd; multiIndId++)
	    {
	      listClusterEl[listClusterElId].indIdCount[multiIndId]=0;
	    }
	  clusterElRead_Single.sizeOfCluster=0;
		
	
	  
	}else{

	//fscanf(fpCluster,"%s %i %i %i %f %i %s %c %c ",chroName1, &startPos, &stopPos, &SVtype, &editProbDist, &editDist, libName, &orient1, &orient2);
	fscanf(fpCluster,"%s %i %s %i %i %g %i %s %c %c ",chroName1, &startPos, chroName2, &stopPos, &SVtype, &editProbDist, &editDist, libName, &orient1, &orient2);
	{//For Mobile Insertions (chrN) the startPos is for reference genome and stopPos is for chrN
	  /////////////////////////////ADD THE CLUSTER TO THE READ////////////////////////
	  readId=findReadName(readNameStr);
	  listReadEl[readId].indId=findIndId(multiLibsCount, libName);
	  listReadEl[readId].libId=findLibId(multiLibsCount, libName);
	  clusterIdElNew = (clusterIdEl *) malloc(sizeof(clusterIdEl));
	  clusterIdElNew->clusterId=listClusterElId;
	  clusterIdElNew->next=listReadEl[readId].next;
	  listReadEl[readId].next=clusterIdElNew;
	  listReadEl[readId].indId=findIndId(multiLibsCount, libName);
	  //////////////////////////////ADD THE READ TO THE CLUSTER/////////////////////////
	  listClusterEl[listClusterElId].clusterId=listClusterElId;
	  clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].readId=readId;
	  clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].indId=findIndId(multiLibsCount, libName);
	  clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].editDistance=editDist;
	  if (weightedHeuristicFlag==1)
	    clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].probEditDist=editProbDist;
	  else 
	    clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].probEditDist=1; // The unweighted verison this score for every mapping is equal to 1
	  clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].posMapLeft=startPos;
	  clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].posMapRight= stopPos;
	  clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].orient1=orient1;
	  clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].orient2=orient2;
	  strcpy(clusterElRead_Single.readMappingElArray[clusterElRead_Single.sizeOfCluster].chroName, chroName1);
	  clusterElRead_Single.sizeOfCluster++;
	}
      }
    }

  sizeListClusterEl=listClusterElId;

}

float scoreForEditDistanceFunc(int totalEditDist)
{
  return totalEditDist*mismatchPenalty;
}


float calWeight(int clusterId, int *countBestSetPicked)
{

  if (listClusterEl[clusterId].oldBestIsGood==1)
    {		
      int count=0;
      for ( count=0; count<totalNumInd; count++)
	countBestSetPicked[count] = listClusterEl[clusterId].bestReadToRemove[count];
      return listClusterEl[clusterId].oldBestScore;
    }
  int idIndPickedAlready=0; // the id of set of individuals which have been selected in previous rounds to have this SV 
  int idIndCanBePicked=0; // the id of individuals which we are considering to have this SV in this round
  int countBestSubset=0; // the number of top set of paired-end reads which had best score till that iteration
  int editDistanceSum[totalNumInd]; // the total sume of edit distance for paired-end reads cconsiderd till now

  //int editDistanceSum; //total sume of edit distance for paired-end reads 

  int supOfIndSeen[totalNumInd]; //the number of paired-end read supports for each individual seen in this round
  double weightedSupOfIndSeen[totalNumInd];// a normalized numer of paired-end read supports for each individual seen 
  //int totalSupSeen; // the number of paired-end read supporting the SV
  int totalEditDist=0; // total of edit distance of all paired-end reads for individuals which we are consiering
  float bestScore, weightNew, normalizedSup=0, normalizedWeightedSup=0;
  int numReadCanBeCovered=0;

  int indCount=0;

  for ( indCount=0; indCount<totalNumInd; indCount++)
    {
      supOfIndSeen[indCount]=0;
      editDistanceSum[indCount]=0;
      countBestSetPicked[indCount]=0;
      weightedSupOfIndSeen[indCount]=0;
	
    }


  bestScore=inf;

  readMappingEl *ptrReadMapping;
  ptrReadMapping=listClusterEl[clusterId].next;

  //printf("%i %i %f\n", clusterId, idIndPickedAlready, weightsForCombination[idIndPickedAlready]);

  while(ptrReadMapping!=NULL) 
    {
      if (listReadEl[ptrReadMapping->readId].readCovered==0 && listClusterEl[clusterId].indIdCount[ptrReadMapping->indId]>-1)
	{
	  supOfIndSeen[ptrReadMapping->indId]++;

	  //////////////////////////////////////THE HUERISTIC//////////////////////////
	  //////////////////////////////////////Instead of taking the total support as number of reads in each cluster, we use summation of 
	  //////////////////////////////////////normalized number of reads. i.e each read contributes as 1/(total number of mappings) it has.
	  ////////////////////////////////////////////////////////////////////////////

	  weightedSupOfIndSeen[ptrReadMapping->indId]=weightedSupOfIndSeen[ptrReadMapping->indId]+ptrReadMapping->probEditDist;
	  if (supOfIndSeen[ptrReadMapping->indId]==minimumSupNeeded)
	    {
	      idIndCanBePicked=idIndCanBePicked+(int)pow(2, ptrReadMapping->indId);			
	      numReadCanBeCovered=numReadCanBeCovered+supOfIndSeen[ptrReadMapping->indId];
	      totalEditDist=totalEditDist+editDistanceSum[ptrReadMapping->indId];
	    } else if (supOfIndSeen[ptrReadMapping->indId]>minimumSupNeeded)
	    {
	      numReadCanBeCovered++;
	      totalEditDist=totalEditDist+ptrReadMapping->editDistance;
	    }
	}

      ptrReadMapping=ptrReadMapping->next;
    }

  weightNew = weightsForCombination[idIndCanBePicked]+scoreForEditDistanceFunc(totalEditDist);
  if (numReadCanBeCovered>0)
    {		
      int count=0;
      for (count=0; count<multiIndCount; count++)
	{
	  if (supOfIndSeen[count]>=minimumSupNeeded)
	    {
	      if (multiIndCount>1)
		{
		  normalizedSup = normalizedSup+((float)(multiIndCount)/(float)(multiIndCount-1))*((float)(totalCoverage-coverageForIndividual[count])/(float)(totalCoverage))*supOfIndSeen[count];
		  normalizedWeightedSup = normalizedWeightedSup + ((float)(multiIndCount)/(float)(multiIndCount-1))*((float)(totalCoverage-coverageForIndividual[count])/(float)(totalCoverage)) * weightedSupOfIndSeen[count];
		} else if (multiIndCount==1)
		{
		  normalizedSup = normalizedSup + supOfIndSeen[count];
		  normalizedWeightedSup = normalizedWeightedSup +  weightedSupOfIndSeen[count];
		}
 
	    }	
	}
      bestScore = ((float)weightNew/(float)normalizedWeightedSup);
    }
  else 
    {
      bestScore = inf;
    }
  int count=0;

  for (count=0; count<multiIndCount; count++)
    {
      if (supOfIndSeen[count]>=minimumSupNeeded)
	{
	  countBestSetPicked[count]=supOfIndSeen[count];
	}	
    }
  listClusterEl[clusterId].oldBestIsGood=1;
  listClusterEl[clusterId].oldBestScore=bestScore;
  for (count=0; count<totalNumInd; count++)
    {
      listClusterEl[clusterId].bestReadToRemove[count]=countBestSetPicked[count];
    }

  return bestScore;
}



int removeInd(int ind, int j)
{
  listClusterEl[j].indIdCount[ind]=-1;
  listClusterEl[j].oldBestIsGood=0;
}


int conflictsBetweenTwoSV_Cord(int startCoord1, int stopCoord1,char SVtype1, int startCoord2, int stopCoord2, char SVtype2)
{
  if ((startCoord1 > stopCoord2 || stopCoord1 < startCoord2))
    return 0;
  if (SVtype1=='A' || SVtype1=='B' || SVtype2=='A' || SVtype2=='B')
    return 0;
  if (SVtype1=='V' && SVtype2=='V')
    {
      if (((startCoord1 > startCoord2) && (stopCoord1 < stopCoord2))||((startCoord1 < startCoord2) && (stopCoord1 > stopCoord2)))
	return 0;
    } else if (SVtype1=='V' && SVtype2=='D')
    {
      if (startCoord1<startCoord2 && stopCoord1>stopCoord2)
	return 0;

    }else if (SVtype1=='D' && SVtype2=='V')
    {
      if (startCoord2<startCoord1 && stopCoord2>stopCoord1)
	return 0;
    }
  return 1;
}

int conflictsAny(int i, int *supInd) // return the individual that SV_i is in conflict with any of the selected SVs
{
  SV_selected *ptrSV;	
  int conflict=0;
  int count=0;
  for (count=0; count<numSV; count++)
    {
      if (strcmp(listSelectedSV[count].chroName, listClusterEl[i].chroName)==0)
	{
	  if(conflictsBetweenTwoSV_Cord(listSelectedSV[count].posStart_SV, listSelectedSV[count].posEnd_SV, listSelectedSV[count].SVtype, listClusterEl[i].posStartSV, listClusterEl[i].posEndSV, listClusterEl[i].SVtype))
	    {
	      ptrSV=listSelectedSV[count].conflict_Next;
	      while (ptrSV!=NULL)
		{
		  if (conflictsBetweenTwoSV_Cord(ptrSV->posStart_SV, ptrSV->posEnd_SV, ptrSV->SVtype, listClusterEl[i].posStartSV, listClusterEl[i].posEndSV, listClusterEl[i].SVtype))
		    {
		      int countInd=0;
		      for (countInd=0; countInd<multiIndCount; countInd++)
			{	
			  if (ptrSV->sup[countInd] > 0 && listClusterEl[i].indIdCount[countInd] <= 0 && listSelectedSV[count].sup[countInd] > 0 && supInd[countInd]>0)
			    {
								
			      removeInd(countInd, i);
			      conflict=1;
			    }
			}
		    }
		  ptrSV=ptrSV->conflict_Next;
		}
	    }
	}
    }
  if (conflict==1) return 1;
  else return -1;
}

int checkConflictNewSelected(SV_selected *ptr1, SV_selected *ptr2, int NotSelectedClus)
{
  int stillCanBeSelected=0;
  int count=0;

  for (count=0; count<multiIndCount; count++)
    {
      if (listClusterEl[NotSelectedClus].indIdCount[count]==0)
	stillCanBeSelected=1;
    }
  if (stillCanBeSelected==0)
    return 0;

  if (conflictsBetweenTwoSV_Cord(ptr2->posStart_SV, ptr2->posEnd_SV, ptr2->SVtype,  listClusterEl[NotSelectedClus].posStartSV, listClusterEl[NotSelectedClus].posEndSV, listClusterEl[NotSelectedClus].SVtype))
    {
      int countInd;
      for (countInd=0; countInd<multiIndCount; countInd++)
	{
	  if (listClusterEl[NotSelectedClus].indIdCount[countInd]==0 && ptr1->sup[countInd]>0 && ptr2->sup[countInd]>0)
	    {
	      listClusterEl[NotSelectedClus].indIdCount[countInd]=-1;
	      listClusterEl[NotSelectedClus].oldBestIsGood=0;
	    }
				
	}
    }
  return 0;
}

int conflictsTwoWay(int i, int j) //returns the individual that SV_i and SV_j are in conflict among the slected SVs
{
  if (strcmp(listSelectedSV[i].chroName, listSelectedSV[j].chroName)!=0)
    return -1;
  if (conflictsBetweenTwoSV_Cord(listSelectedSV[i].posStart_SV, listSelectedSV[i].posEnd_SV, listSelectedSV[i].SVtype, listSelectedSV[j].posStart_SV, listSelectedSV[j].posEnd_SV, listSelectedSV[j].SVtype)==0)
    return -1;
  int count=0;
  for ( count=0; count<multiIndCount; count++)
    {
      if (listSelectedSV[i].sup[count]>0 && listSelectedSV[j].sup[count]>0)
	return count;
    } 
  return -1;
}

int wasNotInConliftNowIsConflict(int i, int j, int *countReads)
// add do the two SV selected i and j not in conflict priousvly 
//however with the additional supported given to cluster "j" they will be in conflict
{

  if (strcmp(listSelectedSV[i].chroName, listSelectedSV[j].chroName)!=0)
    return -1;
  if (conflictsBetweenTwoSV_Cord(listSelectedSV[i].posStart_SV, listSelectedSV[i].posEnd_SV, listSelectedSV[i].SVtype, listSelectedSV[j].posStart_SV, listSelectedSV[j].posEnd_SV, listSelectedSV[j].SVtype)==0)
    return -1;
  int count=0;
  for (count=0; count<multiIndCount; count++)
    {
      if (listSelectedSV[i].sup[count]>0 && listSelectedSV[j].sup[count]>0)
	return -1;
    } 
  for (count=0; count<multiIndCount; count++)
    {
      if (listSelectedSV[i].sup[count]>0 && listSelectedSV[j].sup[count]==0 && countReads[count]>0)
	return count;
    }
  return -1;
}

int addToListOfConflicts(int i, int j, int *countReads) // adds SV i to SV j's conflict list and SV j to SV i's conflict list
{

  SV_selected *newSV;
  newSV=(SV_selected *) malloc(sizeof(SV_selected));
  newSV->posStart_SV=listSelectedSV[i].posStart_SV;
  newSV->posEnd_SV=listSelectedSV[i].posEnd_SV;
  newSV->clusterId=listSelectedSV[i].clusterId;
  newSV->SVtype=listSelectedSV[i].SVtype;
  strcpy(newSV->chroName, listSelectedSV[i].chroName);
  int count=0;
  for ( count=0; count<multiIndCount; count++)
    {
      newSV->sup[count]=listSelectedSV[i].sup[count];
    }		
  newSV->conflict_Next=listSelectedSV[j].conflict_Next;
  listSelectedSV[j].conflict_Next=newSV;
  newSV=listSelectedSV[j].conflict_Next;

}


void addToConflict(int maxWeightSet, int *countReads)// adds the SV numSV to the conflicts tables of previous SVs
{
  int idSV2Add=-1;
  SV_selected *ptrTmp;
  int count=0;

  for ( count=0; count<numSV; count++)
    {
      if (listSelectedSV[count].clusterId==maxWeightSet)
	{
	  idSV2Add=count;
	}

    }
  if (idSV2Add==-1)
    {
      strcpy(listSelectedSV[numSV].chroName, listClusterEl[maxWeightSet].chroName);
      listSelectedSV[numSV].posStart_SV=listClusterEl[maxWeightSet].posStartSV;
      listSelectedSV[numSV].clusterId=maxWeightSet;
      listSelectedSV[numSV].posEnd_SV=listClusterEl[maxWeightSet].posEndSV;
      listSelectedSV[numSV].conflict_Next=NULL;
      listSelectedSV[numSV].SVtype=listClusterEl[maxWeightSet].SVtype;
      int count=0;
      for ( count=0; count<multiIndCount; count++)
	{
	  listSelectedSV[numSV].sup[count]=listClusterEl[maxWeightSet].indIdCount[count];
	}
      int i=0;
      for ( i=0; i<numSV; i++)
	{
	  if (conflictsTwoWay(numSV, i)>=0)
	    {
	      addToListOfConflicts(i, numSV, countReads);
	      addToListOfConflicts(numSV, i, countReads);
	    }	
	}


      for ( count=0; count<sizeListClusterEl; count++)
	{
	  ptrTmp=listSelectedSV[numSV].conflict_Next;
	  if (conflictsBetweenTwoSV_Cord(listSelectedSV[numSV].posStart_SV, listSelectedSV[numSV].posEnd_SV, listSelectedSV[numSV].SVtype, listClusterEl[count].posStartSV, listClusterEl[count].posEndSV, listClusterEl[count].SVtype) && (strcmp(listSelectedSV[numSV].chroName, listClusterEl[count].chroName)==0))
	    {
	      while(ptrTmp!=NULL)
		{	
		  checkConflictNewSelected(&listSelectedSV[numSV], ptrTmp, count);
		  ptrTmp=ptrTmp->conflict_Next;
		}
	    }
	}
      numSV++;
    }else{
    int count=0;
    for (count=0; count<numSV; count++)
      {
	if (count!=idSV2Add && wasNotInConliftNowIsConflict(count, idSV2Add, countReads)>0)
	  {
	    addToListOfConflicts(count, idSV2Add, countReads);
	    addToListOfConflicts(idSV2Add, count, countReads);
	  }
      }
    for ( count=0; count<multiIndCount; count++)
      {
	listSelectedSV[idSV2Add].sup[count]=listSelectedSV[idSV2Add].sup[count]+countReads[count];
      }

    for ( count=0; count<sizeListClusterEl; count++)
      {
	ptrTmp=listSelectedSV[idSV2Add].conflict_Next;
	if(conflictsBetweenTwoSV_Cord(listSelectedSV[idSV2Add].posStart_SV, listSelectedSV[idSV2Add].posEnd_SV, listSelectedSV[idSV2Add].SVtype, listClusterEl[count].posStartSV, listClusterEl[count].posEndSV, listClusterEl[count].SVtype) && (strcmp(listSelectedSV[idSV2Add].chroName, listClusterEl[count].chroName)==0))
	  {
	    while(ptrTmp!=NULL)
	      {
		checkConflictNewSelected(&listSelectedSV[numSV], ptrTmp, count);
		ptrTmp=ptrTmp->conflict_Next;
	      }
	  }
      }

  }


}


int search_MobileName(int pos)
{
  int count=0;
  for( count=0; count<listMobileElSize; count++)
    {
      if (pos<= listMobileEl[count].endPos+10 && pos >= listMobileEl[count].startPos-10)
	return (count);
    }
}


float avgEditDistReadSupportingCluster(readMappingEl *list)
{
  int countTotal=0, editDistanceSum=0;
	
  while(list!=NULL)
    {
      editDistanceSum = editDistanceSum + list->editDistance;
      countTotal++;
      list=list->next;
    }

  return ((float)(editDistanceSum)/(float)countTotal);

}



float avgEditDistIndReadSupportCluster(readMappingEl *list, int countInd)
{
  int countTotal=0, editDistanceSum=0;
  while(list!=NULL)
    {
      if (list->indId==countInd)
	{
	  editDistanceSum = editDistanceSum + list->editDistance;
	  countTotal++;
	}
      list=list->next;
    }
  if (countTotal==0)
    return 0; else
    return ((float)(editDistanceSum)/(float)countTotal);
}



int outputReadsSupportingCluster(readMappingEl *list, FILE *fpOut)
{
  while(list!=NULL)
    {
      //		fprintf(fpOut, "%s %s %i %i %i %g %c %c %s %s\n",listReadEl[list->readId].readName, list->chroName, list->posMapLeft, list->posMapRight, list->editDistance, list->probEditDist, list->orient1, list->orient2, multiLibs[listReadEl[list->readId].libId].libName, multiInd[listReadEl[list->readId].indId]);
      list=list->next;
    }
}

double * calculateTheHeuristicScore(readMappingEl *list)
{
  double *result;
  result = (double *) malloc(multiIndCount*sizeof(double));
  int count=0;
  for( count=0; count<multiIndCount; count++)
    {
      result[count]=0;
    }	

  while(list!=NULL)
    {
      result[listReadEl[list->readId].indId] = result[listReadEl[list->readId].indId] + list->probEditDist;
      list=list->next;
    }
  return result;

}

int outputCluster(int set, FILE *fpOut)
{
  readMappingEl *ptrReadMappingSelected;
  int outInsLeft=maxChroSize, inInsLeft=0, inInsRight=maxChroSize, outInsRight=0;
  char mobileName[strSize];
  int idMobile;
  int supTotal=0;
  int F_Side=0;
  int R_Side=0;
  double *heuristicScore;
  int count=0;
  for ( count=0; count<multiIndCount; count++)
    {
      if (listClusterEl[set].indIdCount[count]>0)
	supTotal=supTotal+listClusterEl[set].indIdCount[count];
    }
  if (listClusterEl[set].SVtype=='A')
    {
      ptrReadMappingSelected = listClusterEl[set].readMappingSelected;
      while(ptrReadMappingSelected!=NULL)
	{
	  if (ptrReadMappingSelected->orient1=='F' && ptrReadMappingSelected->orient2=='R')
	    {
	      if (ptrReadMappingSelected->posMapRight - multiLibs[listReadEl[ptrReadMappingSelected->readId].libId].readLen < outInsLeft)
		outInsLeft = ptrReadMappingSelected->posMapRight - multiLibs[listReadEl[ptrReadMappingSelected->readId].libId].readLen;
	      if (ptrReadMappingSelected->posMapRight -  multiLibs[listReadEl[ptrReadMappingSelected->readId].libId].readLen > inInsLeft)
		inInsLeft = ptrReadMappingSelected->posMapRight - multiLibs[listReadEl[ptrReadMappingSelected->readId].libId].readLen;
	      ;
	      F_Side=1;
	    }else if (ptrReadMappingSelected->orient1=='R' && ptrReadMappingSelected->orient2=='F')
	    {
	      if (ptrReadMappingSelected->posMapRight < inInsRight)
		inInsRight = ptrReadMappingSelected->posMapRight;
	      if (ptrReadMappingSelected->posMapRight > outInsRight)
		outInsRight = ptrReadMappingSelected->posMapRight;
	      R_Side=1;
	    }
	  ptrReadMappingSelected=ptrReadMappingSelected->next;
	}
      strcpy(mobileName, listMobileEl[search_MobileName(listClusterEl[set].next->posMapRight)].nameMobile);
      //		printf("%s %i\n", mobileName, search_MobileName(listClusterEl[set].next->posMapRight));
      idMobile = search_MobileName(listClusterEl[set].next->posMapRight);
      outInsLeft = outInsLeft - listMobileEl[idMobile].startPos;
      inInsLeft = inInsLeft - listMobileEl[idMobile].startPos;
      inInsRight = inInsRight - listMobileEl[idMobile].startPos;
      outInsRight = outInsRight - listMobileEl[idMobile].startPos; 


      if (F_Side == 0)
	{
	  outInsLeft=0;
	  inInsLeft=0;
	}else if (R_Side==0)
	{	inInsRight = listMobileEl[idMobile].endPos - listMobileEl[idMobile].startPos;
	  outInsRight = listMobileEl[idMobile].endPos - listMobileEl[idMobile].startPos;
			
	}
      if (outInsLeft<0)
	outInsLeft=0;
      if (inInsLeft<0)
	inInsLeft=0;
      if (inInsRight<0)
	inInsRight=0;
      if (outInsRight<0)
	outInsRight=0;

    }else if (listClusterEl[set].SVtype=='B')
    {
		
      ptrReadMappingSelected = listClusterEl[set].readMappingSelected;
      while(ptrReadMappingSelected!=NULL)
	{
	  if (ptrReadMappingSelected->orient1=='F' && ptrReadMappingSelected->orient2=='F')
	    {
	      if (ptrReadMappingSelected->posMapRight > outInsRight)
		outInsRight = ptrReadMappingSelected->posMapRight;
	      if (ptrReadMappingSelected->posMapRight < inInsRight)
		inInsRight = ptrReadMappingSelected->posMapRight;
	      F_Side = 1;
	    }else if (ptrReadMappingSelected->orient1=='R' && ptrReadMappingSelected->orient2=='R')
	    {
	      if (ptrReadMappingSelected->posMapRight - multiLibs[listReadEl[ptrReadMappingSelected->readId].libId].readLen > inInsLeft)
		inInsLeft = ptrReadMappingSelected->posMapRight - multiLibs[listReadEl[ptrReadMappingSelected->readId].libId].readLen;
	      if (ptrReadMappingSelected->posMapRight- multiLibs[listReadEl[ptrReadMappingSelected->readId].libId].readLen < outInsLeft)
		outInsLeft = ptrReadMappingSelected->posMapRight - multiLibs[listReadEl[ptrReadMappingSelected->readId].libId].readLen;
	      R_Side = 1;
	    }
	  ptrReadMappingSelected=ptrReadMappingSelected->next;
	}
      strcpy(mobileName, listMobileEl[search_MobileName(listClusterEl[set].next->posMapRight)].nameMobile);
      //printf("%s %i\n", mobileName, search_MobileName(listClusterEl[set].next->posMapRight)); 
      idMobile = search_MobileName(listClusterEl[set].next->posMapRight);
      outInsLeft = outInsLeft - listMobileEl[idMobile].startPos;
      inInsLeft = inInsLeft - listMobileEl[idMobile].startPos;
      inInsRight = inInsRight - listMobileEl[idMobile].startPos;
      outInsRight = outInsRight - listMobileEl[idMobile].startPos; 
      if (R_Side == 0)
	{
	  outInsLeft=0;
	  inInsLeft=0;
	}else if (F_Side==0)
	{	inInsRight = listMobileEl[idMobile].endPos - listMobileEl[idMobile].startPos;
	  outInsRight = listMobileEl[idMobile].endPos - listMobileEl[idMobile].startPos;
			
	}
      if (outInsLeft<0)
	outInsLeft=0;
      if (inInsLeft<0)
	inInsLeft=0;
      if (inInsRight<0)
	inInsRight=0;
      if (outInsRight<0)
	outInsRight=0;
    }



  heuristicScore = calculateTheHeuristicScore(listClusterEl[set].readMappingSelected);


  if (listClusterEl[set].SVtype=='A')
    {
	

      fprintf(fpOut,"Chr:%s Start_Outer:%i Start_Inner:%i End_Inner:%i End_Outer:%i SVtype:%c MobileName:%s StartLeftIns:%i EndLeftIns:%i StartRightIns:%i EndRightIns:%i sup:%i Sum_Weight:0 AvgEditDist:%f", listClusterEl[set].chroName, listClusterEl[set].posStartSV_Outer, listClusterEl[set].posStartSV, listClusterEl[set].posEndSV, listClusterEl[set].posEndSV_Outer, listClusterEl[set].SVtype, mobileName,outInsLeft, inInsLeft, inInsRight, outInsRight, supTotal,avgEditDistReadSupportingCluster( listClusterEl[set].readMappingSelected  ) );
    } else if (listClusterEl[set].SVtype=='B')
    {
      fprintf(fpOut,"Chr:%s Start_Outer:%i Start_Inner:%i End_Inner:%i End_Outer:%i SVtype:%c MobileName:%s StartLeftIns:%i EndLeftIns:%i StartRightIns:%i EndRightIns:%i sup:%i Sum_Weight:0 AvgEditDist:%f", listClusterEl[set].chroName, listClusterEl[set].posStartSV_Outer, listClusterEl[set].posStartSV, listClusterEl[set].posEndSV, listClusterEl[set].posEndSV_Outer, listClusterEl[set].SVtype, mobileName,outInsLeft, inInsLeft, inInsRight, outInsRight, supTotal,avgEditDistReadSupportingCluster( listClusterEl[set].readMappingSelected  ) );
    }
  else
    {
      fprintf(fpOut,"Chr:%s Start_Outer:%i Start_Inner:%i End_Inner:%i End_Outer:%i SVtype:%c sup:%i Sum_Weight:0 AvgEditDits:%f", listClusterEl[set].chroName, listClusterEl[set].posStartSV_Outer, listClusterEl[set].posStartSV, listClusterEl[set].posEndSV, listClusterEl[set].posEndSV_Outer, listClusterEl[set].SVtype, supTotal, avgEditDistReadSupportingCluster( listClusterEl[set].readMappingSelected  ));
    }
  int countInd;
  for (countInd=0; countInd<multiIndCount; countInd++)
    {
      fprintf(fpOut, " Lib:%s LibSup:%i LibHurScore:%g AvgEditDistInd:%g", multiInd[countInd], listClusterEl[set].indIdCount[countInd], heuristicScore[countInd], avgEditDistIndReadSupportCluster(listClusterEl[set].readMappingSelected, countInd));
    }

  if (listClusterEl[set].SVtype=='D')
    {
      fprintf(fpOut, " minDelLen:%i maxDelLen:%i", listClusterEl[set].minDelLength, listClusterEl[set].maxDelLength);
    }

  fprintf(fpOut, "\n");

  outputReadsSupportingCluster(listClusterEl[set].readMappingSelected, fpOut);
  free(heuristicScore);

}

int markReadsCovered(int clusterId, int *countReads)
{
  int minDelLength=1000000000, maxDelLength=1000000000; //Only will be used for the type deletion

  readMappingEl *ptrRead;
  clusterIdEl * ptrClusterId;
  ptrRead=listClusterEl[clusterId].next;
  int totalReadRemove=0;
  int readsToMark[totalNumInd];
  int count=0;
  for ( count=0; count<multiIndCount; count++)
    {
      totalReadRemove=totalReadRemove+countReads[count];
      readsToMark[count]=countReads[count];
    }

  //This part will mark the reads covered and calculate the minimum and maximum length of deletion 
	
	

  while(ptrRead!=NULL && totalReadRemove>0)
    {
      if (readsToMark[ptrRead->indId]>0 && listReadEl[ptrRead->readId].readCovered==0)
	{
	  ptrClusterId=listReadEl[ptrRead->readId].next;
	  listClusterEl[clusterId].indIdCount[ptrRead->indId]++;
			
			
	  listClusterEl[clusterId].readMappingSelected=addToListOfReads(listClusterEl[clusterId].readMappingSelected, *ptrRead); 
			
	  if (listClusterEl[clusterId].SVtype=='D' && minDelLength > ptrRead->posMapRight - ptrRead->posMapLeft - multiLibs[listReadEl[ptrRead->readId].libId].maxInstSize)
	    {
	      minDelLength=ptrRead->posMapRight - ptrRead->posMapLeft - multiLibs[listReadEl[ptrRead->readId].libId].maxInstSize;
			
	    }
			
	  if (listClusterEl[clusterId].SVtype=='D' && maxDelLength > ptrRead->posMapRight - ptrRead->posMapLeft - multiLibs[listReadEl[ptrRead->readId].libId].minInstSize)
	    {
	      maxDelLength=ptrRead->posMapRight - ptrRead->posMapLeft - multiLibs[listReadEl[ptrRead->readId].libId].minInstSize;
			
	    }
	
	  listReadEl[ptrRead->readId].readCovered=1;
	  readsToMark[ptrRead->indId]--;
	  while (ptrClusterId!=NULL)
	    {

	      listClusterEl[ptrClusterId->clusterId].oldBestIsGood=0;
	      ptrClusterId=ptrClusterId->next;
	    }
	  totalReadRemove--;
	}
      ptrRead=ptrRead->next;
    }

  if (listClusterEl[clusterId].SVtype=='D')
    {
      listClusterEl[clusterId].minDelLength=minDelLength;
      listClusterEl[clusterId].maxDelLength=maxDelLength;
    }
	

}

int bufferIsUseful()
{

  float bestWeight=inf;
  int bestWeightId=-1;

  int count=0;
  for ( count=0; count<countInBuffer; count++)
    {
      if (listClusterEl[listClusterInBuffer[count].clusterId].oldBestIsGood==0)
	{
	  listClusterEl[listClusterInBuffer[count].clusterId].oldBestScore=calWeight(listClusterInBuffer[count].clusterId, listClusterEl[listClusterInBuffer[count].clusterId].bestReadToRemove );
	  listClusterEl[listClusterInBuffer[count].clusterId].oldBestIsGood=1;
	  listClusterInBuffer[count].score=listClusterEl[listClusterInBuffer[count].clusterId].oldBestScore;
	}
     
    }
 
  for ( count=0; count<countInBuffer; count++)
    {
      if (bestWeight > listClusterEl[listClusterInBuffer[count].clusterId].oldBestScore)
	{
	  bestWeight=listClusterEl[listClusterInBuffer[count].clusterId].oldBestScore;
	  bestWeightId=count;	
	}
     
    }
 
  if(bestWeight<=maxScoreInBuffer && bestWeight!=inf)
    return 1;
  else return 0;
}


int pickSet()
{
  float bestWeight=inf;
  float newWeight;
  int countReads[totalNumInd];
  int bestReads[totalNumInd];
  int bestWeightSet=0;
  int setToCheckId=0;
  int bestReadCovered;
  maxScoreInBuffer=inf-1;
  while(numCallsRequested>0)
    {
      bestWeight=inf;
      bestWeightSet=-1;
     
      if(!bufferIsUseful())
	{
	  emptyBuffer();
	  int clusterCounter=0;
	  for (clusterCounter=0; clusterCounter<sizeListClusterEl; clusterCounter++)
	    {	
	      newWeight=calWeight(clusterCounter, countReads);
	      addToBuffer(newWeight, clusterCounter);
	      if (bestWeight > newWeight)
		{	
		  bestWeightSet = clusterCounter;
		  int count=0;
		  for ( count=0; count<multiIndCount; count++)
		    {
		      bestReads[count]=countReads[count];
		    }
		  bestWeight=newWeight;
		}
	    }

			
	}else{
	bestWeightSet=bestFromBuffer();
	bestWeight=listClusterEl[bestWeightSet].oldBestScore;
	int count=0;
	for ( count=0; count<multiIndCount; count++)
	  {
	    bestReads[count]=listClusterEl[bestWeightSet].bestReadToRemove[count];
	  }
      }
	
      if (bestWeight<inf)
	{
	  if (conflictResFlag == 0 || conflictsAny(bestWeightSet, bestReads)==-1)
	    {


	      {
		markReadsCovered(bestWeightSet, bestReads);
	
		if (conflictResFlag==1)
		  addToConflict(bestWeightSet, bestReads);
	
		numCallsRequested--;
		listClusterEl[bestWeightSet].oldBestIsGood=0;
	      }
	    }
	}
      else 
	{
	  /*
	    printf("aaaaaaa%f\n", bestWeight);
	  */
	  return 0;
	}
		
    }
}




int outputPickedCluster(FILE *fpOut)
{

  int totalSup=0;
  int count=0;
  for ( count=0; count<sizeListClusterEl; count++)
    {
      totalSup=0;
      int count2=0;
      for (count2=0; count2<multiIndCount; count2++)
	{
	  totalSup=totalSup+listClusterEl[count].indIdCount[count2];
	}
      if (totalSup>0)
	{
	  outputCluster(count, fpOut);
	}
			
    }
}



void readMobileElements(FILE *fp)
{
  char nameM[50], nameC[50];
  int pos1, pos2; 
  
  while (fscanf(fp,"%s\t%i\t%i\t%s\n", nameC,&pos1, &pos2, nameM)!=EOF)
    {
      strcpy(listMobileEl[listMobileElSize].nameMobile, nameM);
      listMobileEl[listMobileElSize].startPos=pos1;
      listMobileEl[listMobileElSize].endPos=pos2;
      listMobileElSize++;
    }
}


int vh_setcover(bam_info **in_bams, int num_bams, char* outputread, char* outputfile, char* svfile){
  FILE *readFp=NULL, *clusterFp=NULL, *libFp=NULL, *weightsFp=NULL, *fpOut=NULL, *fpMobile=NULL, *coverageFp=NULL;
  numCallsRequested = 10000; // have to fix

  readFp=safe_fopen(outputread,"r");
  clusterFp=safe_fopen(outputfile,"r");
  fpOut=safe_fopen(svfile,"w");

  init(in_bams, num_bams, readFp, clusterFp, weightsFp, coverageFp);
  if (fpMobile!=NULL)
    readMobileElements(fpMobile);
  pickSet();	
  outputPickedCluster(fpOut);
}

