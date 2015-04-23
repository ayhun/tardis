#include <string.h>
#include "vh_main.h"
#include "vh_divethandler.h"
#include "vh_maximalCluster.h"
#include "vh_gaphandler.h"
#include "vh_repeathandler.h"

FILE *fileOutput = NULL;
char **allReadNameList;

void vh_quitProgram (int exitCode)
{
  if (exitCode != 0)
    {
      if (!strcmp (g_error_message, ""))
	strcpy (g_error_message, "Unexpected Error Found!\n");

      char errMsg[150];
      sprintf (errMsg, "Program Exited with Errors:%d", exitCode);
      vh_logError (errMsg);
      switch (exitCode)
	{
	case EXIT_CODE_ARG_ERROR:
	  sprintf (errMsg, "Error in options: %s", g_error_message);
	  vh_logError (errMsg);
	  break;
	case EXIT_CODE_DIVET_ERROR:
	  sprintf (errMsg, "Error in DIVET file: %s", g_error_message);
	  vh_logError (errMsg);
	  break;
	case EXIT_CODE_MEMORY_ERROR:
	  sprintf (errMsg, "Memory Problem Occured: %s", g_error_message);
	  vh_logError (errMsg);
	  break;
	default:
	  sprintf (errMsg, "Uncategorized Error Found: %s", g_error_message);
	  vh_logError (errMsg);
	  break;
	}
      exit (exitCode);
    }
  else
    {
      vh_logInfo ("Program Exited Successfully.\r\n");
      exit (0);
    }
}


void vh_pruneAndNormalizeDivets (struct LibraryInfo *lib, double preProsPrune,
			 int overMapLimit)
{
  struct DivetRow *cursor = lib->head->next;



  if (cursor == NULL)
    return;
  //  printf ("%f\n", preProsPrune);


  //We are not checking the first and last one
  while (cursor != NULL && cursor->next != NULL)
    {
      if (((cursor->next->phredScore /
	    cursor->next->readName->sumPhredValue) < preProsPrune)
	  || (cursor->next->readName->occurrences > overMapLimit))
	{
	  struct DivetRow *toBeDeleted = cursor->next;
	  cursor->next = cursor->next->next;
	  lib->size--;
	  free (toBeDeleted->chroName);
	  free (toBeDeleted);
	}
      else
	{
	  cursor->next->phredScore =
	    cursor->next->phredScore / cursor->next->readName->sumPhredValue;
	  cursor = cursor->next;
	}


    }
}


static int vh_cmprReadNameStr (const void *a, const void *b)
{
  return strcmp (*(char **) a, *(char **) b);
}



//We assume that the arguments are valid at this point and the file exists
//void vh_clustering (char *libFileAdrs, bam_info* in_bam, char *gapFileName,
//  char *repeatFileName, double preProsPrune, char *outputFile, char *outputRead, int overMapLimit)
//{

void vh_clustering (bam_info* in_bam, char *gapFileName,
  char *repeatFileName, double preProsPrune, char *outputFile, char *outputRead, int overMapLimit)
{

  int totalNumUniqueReads = 0;
  int indexStart = 0;
  int count;
  int i,j;
  struct LibraryInfo *newLibInfo, *cursor, *t;

  vh_readInitFile ();
  vh_readChros (in_bam);
  vh_readGapTable (gapFileName);
  vh_readRepeatTable (repeatFileName);


  g_libInfo = NULL;

  for (i=0;i<in_bam->num_libraries;i++)
    {
      newLibInfo = (struct LibraryInfo *) malloc (sizeof (struct LibraryInfo));
      strcpy (newLibInfo->libName, in_bam->libraries[i]->libname);
      strcpy (newLibInfo->libFileAdrs, in_bam->libraries[i]->divet);
      newLibInfo->minDelta = in_bam->libraries[i]->conc_min - 2 * in_bam->libraries[i]->read_length;
      newLibInfo->maxDelta = in_bam->libraries[i]->conc_max - 2 * in_bam->libraries[i]->read_length;
      if (newLibInfo->minDelta < 0)
	newLibInfo->minDelta = 0;
      if (newLibInfo->maxDelta < 0)
	newLibInfo->maxDelta = 0;
      newLibInfo->readLen = in_bam->libraries[i]->read_length;
      newLibInfo->hash =(struct ReadName **) malloc (NHASH * sizeof (struct ReadName *));
      for (j = 0; j < NHASH; j++)
	newLibInfo->hash[j] = NULL;
      newLibInfo->head = NULL;
      newLibInfo->tail = NULL;
      newLibInfo->size = 0;
      newLibInfo->next = NULL;
      
      //If first one in the linkedlist                                                                                                                                                                          
      if (g_libInfo == NULL)
	{
	  g_libInfo = newLibInfo;
	}
      else                  //add to the end of the linked list                                                                                                                                                 
	{
	  for (t = g_libInfo; t->next != NULL; t = t->next)
	    ;        //Skip till end of LinkedList                                                                                                             	  
	  t->next = newLibInfo;
	}
    }
  

  cursor = g_libInfo;
  fileOutput = fopen (outputFile, "w");
  //printf ("prune %f\n", preProsPrune);

  for (; cursor; cursor = cursor->next)
    {
      vh_logInfo ("Reading Divet Files ...");
      vh_loadDivetFile (cursor);
      sprintf (g_loggerMsgBuffer, "%d rows loaded successfully.",cursor->size);
      vh_logInfo (g_loggerMsgBuffer);
      vh_pruneAndNormalizeDivets (cursor, preProsPrune, overMapLimit);
      sprintf (g_loggerMsgBuffer, "%d rows after pruning.", cursor->size);
      vh_logInfo (g_loggerMsgBuffer);
    }

  if (strcmp (outputRead, ""))
    {
      vh_logInfo ("Writing ReadName Sorted");
      cursor = g_libInfo;
      FILE *fileOutputReadName;
      fileOutputReadName = fopen (outputRead, "w");
      while (cursor != NULL)

	{
	  totalNumUniqueReads =
	    totalNumUniqueReads + vh_countNumReads (cursor->hash);
	  cursor = cursor->next;
	}
      cursor = g_libInfo;
      allReadNameList =	(char **) malloc (totalNumUniqueReads * sizeof (char *));
      for (; cursor; cursor = cursor->next)
	{
	  indexStart = vh_exportToArray(cursor->hash, allReadNameList, indexStart);
	}
      qsort (allReadNameList, totalNumUniqueReads, sizeof (char *),vh_cmprReadNameStr);
      fprintf (fileOutputReadName, "%i\n", totalNumUniqueReads + 1);
      for (count = 0; count < totalNumUniqueReads; count++)
	{
	  fprintf (fileOutputReadName, "%s\n", allReadNameList[count]);
	}
      fclose (fileOutputReadName);
    }


  //TODO: Someone should free the memory allocated for the divets and for the libraryinfos

  //Just for test
  /*LibraryInfo* l = g_libInfo;
     for (; l; l = l->next)
     {
     printf("Lib: %s, Adrs: %s, Min: %d, Max: %d, Len: %d, size: %d\n", l->libName, l->libFileAdrs, l->minDelta, l->maxDelta, l->readLen, l->size);
     DivetRow* tmp = l->head;
     for (int counter = 0;tmp != NULL && counter < 10; tmp = tmp->next, counter++)
     printDivet(tmp);
     printf("--------\n");

     char* nameToSearch = "IL22_348:5:285:83:287";
     ReadName* r = getReadNameFromHash(l->hash, nameToSearch);
     if (r != NULL)
     sprintf(g_loggerMsgBuffer, "%s\t%d\t%f\t%f\n", r->readName, r->occurrences, r->minEdit, r->meanPhredValue);
     else
     sprintf(g_loggerMsgBuffer, "%s Not Found!\n", nameToSearch);
     logOutput(g_loggerMsgBuffer);
     printf("=============\n");
     }
   */
  for (i = 0; i < g_chroTableSize; i++)
    {
      fprintf(stderr, "\r                                                        ");
      fflush(stderr);
      fprintf(stderr, "\rProcessing chromosome %s ", g_chroTable[i].chroName);
      fflush(stderr);
      vh_initializeReadMapping_Deletion (g_chroTable[i].chroName,
				      g_chroTable[i].size);
      fprintf(stderr, ".");
      fflush(stderr);      
      vh_createDeletionClusters (g_chroTable[i].size);
      fprintf(stderr, ".");
      fflush(stderr);      
      vh_finalizeReadMapping (g_chroTable[i].chroName, g_chroTable[i].size);
      fprintf(stderr, ".");
      fflush(stderr);      
      vh_initializeReadMapping_Inversion (g_chroTable[i].chroName,g_chroTable[i].size);
      fprintf(stderr, ".");
      fflush(stderr);      
      vh_createInversionClusters (g_chroTable[i].size);
      fprintf(stderr, ".");
      fflush(stderr);      
      vh_finalizeReadMapping (g_chroTable[i].chroName, g_chroTable[i].size);
      fprintf(stderr, ".");
      fflush(stderr);      
      vh_initializeReadMapping_Insertion (g_chroTable[i].chroName,g_chroTable[i].size);
      fprintf(stderr, ".");
      fflush(stderr);      
      vh_createInsertionClusters (g_chroTable[i].size);
      fprintf(stderr, ".");
      fflush(stderr);      
      vh_finalizeReadMapping (g_chroTable[i].chroName, g_chroTable[i].size);
      fprintf(stderr, ".");
      fflush(stderr);      
    }
  fprintf(stderr, "\n");
  fclose (fileOutput);


  /* free g_libInfo */
  cursor = g_libInfo;
  t = cursor;
  while (t != NULL)
  {
    t = cursor->next;
    /* FEREYDOUN: Free the hash as well. This is likely a linked list */
    for (i = 0; i < NHASH; i++)
      if (cursor->hash != NULL)
	free(cursor->hash[i]);
    free(cursor->hash);
    free(cursor);
  }  


}


int vh_isValid ()
{
  //TODO: To be completed. Check for files to exists, etc.

  return 1;
}


