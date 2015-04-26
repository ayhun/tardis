/*
 *  weighted_SetCover_MultiColor_W_ConflictResghted_SetCover_MultiColor_W_ConflictRes2_W_EditEdit.cpp2.h
 *  
 *
 *  Created by Fereydoun Hormozdiari on 24/10/11.
 *  Copyright 2011 Simon Fraser University. All rights reserved.
 *
 */
#ifndef __WEIGHTED_SETCOVER_MULTICOLOR__
#define __WEIGHTED_SETCOVER_MULTICOLOR__

#include <stdlib.h>
#include <stdio.h>
#include "../processbam.h"
#include "../common.h"




#define  maxChroSize 1000000000
#define  inf 10000000 // a constant to rpersent infinity
#define  maxLengthSV_Del 2000000 // maximum length of SV allowed
#define  maxLengthSV_Inv 20000000
#define  maxListClusterSize 300000 // maximum size of a cluster
#define  maxClustersAllowed 7000000 //maximum number of clusters to read
#define  strSize 100 //maximum size of a constant string, such as readName length
#define  maxNumSV 1000000 //maximum number of SV's allowed to be reported.
#define  totalNumInd 100
//#define  maxTotalNumLibs 200
#define  TRUE 1
#define FALSE 0
#define true 1
#define false 0
 
typedef struct multiLib{
  char libName[strSize];
  int indId; // index to the multiInd[X] 
  int maxInstSize;
  int minInstSize;
  int readLen;
}multiLib;

multiLib *multiLibs; ///[maxTotalNumLibs];
//int multiLibsCount;

typedef struct SV_selected{ // SVs which are selected for output are kept here for conflict resolution 
  char chroName[strSize];
  int clusterId; // ID of the cluster selected
  char SVtype;// D: deletion, V: Inversion, I: insertion
  int posStart_SV, posEnd_SV;// The inside coordinates
  int posStart_SV_Outer, posEnd_SV_Outer;//The outside coordinates
  int sup[totalNumInd]; // support for the SV picked for each individual
  struct SV_selected* conflict_Next; // keep a link list of all other SVs selected which are in conlift with this SV (in hapolid genome)
} SV_selected;


SV_selected listSelectedSV[maxNumSV];// the array of all the SVs selected till now 
int numSV; // total number of distinct SVs picked

////////////////////////////////////////////////THE DATA STRUCTURES TO HOLD INFORMATION FOR EACH MAPPING


typedef struct clusterIdEl{ // used to create alinked list for cluster ids
  int clusterId;
  struct clusterIdEl *next;
}clusterIdEl;




typedef struct readEl{ // holds the information for each read
  char readName[strSize]; 
  int readId;
  int indId; 
  int libId;
  int readCovered;//0 not covered, 1 covered;
  struct clusterIdEl *next; // a link list of the clusters which have this read as support
} readEl;



struct readEl *listReadEl; // Array of all reads 
int sizeListReadEl;// number of total reads



///////////////////////////////////////////////THE DATA STRUCTURES FOR HOLDING INFORMATION FOR EACH CLUSTERS


typedef struct readMappingEl{ //  Link list for the information of each Read Mapping for each cluster
  int readId; //read id [0 to sizeListReadEl] is an unique identifier for each read (and index to array listReadEL).
  char chroName[strSize];
  int posMapLeft;
  int posMapRight;
  int indId;
  int editDistance;
  char orient1;
  char orient2;
  float probEditDist;
  struct readMappingEl *next; //ptr to the next readMappingEl in this cluster 
}readMappingEl;

typedef struct clusterEl{// holds the informtion for each cluster of paired-end reads
  int clusterId;
  char chroName[strSize];
  int posStartSV;
  int posEndSV;
  int posStartSV_Outer;
  int posEndSV_Outer;
  int minDelLength;// Only used for deletion. Repersents the minimum size of deletion predicted
  int maxDelLength;//Only used for deletion. Repersents the maximum size of deletion predicted
  char SVtype;//V: inversion, D: Deletion, I: insertion
  int indIdCount[totalNumInd]; //If this cluster is picked as one of the SVs it shows the number of supporting paired-end reads selected for each individual (0 : means that we have not picked any support for that inidividual for this SV. -1: means that for this individual this SV is in conflict with SVs picked before - NEVER PICK AN SV FOR AN IND WITH SUP -1 ).
  struct readMappingEl *next; //a array of number of individuals for each list of read mappings (i.e. for each individual we have adifferent list of read mappings and they should be sorted).
  struct readMappingEl *readMappingSelected; // The link list of all the read mappings which have been selected by set cover for this SV.
  int oldBestIsGood; // is the last best score computed for this SV still the best or things have changes
  float oldBestScore; // Whats is the best old score
  int bestReadToRemove[totalNumInd]; // the support for each individual for the best old score(if oldBestIsGood is true then this array is still the best selection from this cluster).
  /* indIdCount[totalNumInd] : repersents number of paired-end read support for each individual for this cluster selected and covered till now
     bestReadToRemove[totalNumInd]: if an individual from this cluster wants to be covered or additional set of reads want to be covered this is the best possible selection
  */

}clusterEl;

typedef struct clusterElRead{// the cluster which is being read and processed before puting all inside a large array (temp cluster)
  int clusterId;
  char chroName[strSize];
  char SVtype;
  int sizeOfCluster;
  struct readMappingEl readMappingElArray[maxListClusterSize];
}clusterElRead;

int vh_setcover(bam_info **in_bam, int num_bams, char* outputread, char* outputfile, char* svfile);


#endif
