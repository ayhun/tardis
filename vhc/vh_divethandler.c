#include <time.h>
#include "vh_divethandler.h"
#include "vh_gaphandler.h"
#include "vh_repeathandler.h"

struct LibraryInfo *g_libInfo = NULL;

struct DivetRow *createDivetRow (struct ReadName *hash[], char *readName, char *chroName, char *locMapLeftStart, char *locMapLeftEnd, char *orientationLeft, char *locMapRightStart, char *locMapRightEnd, char *orientationRight, char *svType, char *editDistance, char *avgQual,	//skip
				 char *phredScore, struct LibraryInfo *libInfo, int id)
{
  struct DivetRow *newRow = (DivetRow *) malloc (sizeof (DivetRow));
  if (newRow == NULL)
    {
      strcpy (g_error_message, "Memory Problem in loading divet file!\n");
      vh_quitProgram (EXIT_CODE_MEMORY_ERROR);
    }

  newRow->readName = NULL;
  newRow->next = NULL;

  int len;

  len = strlen (chroName);
  newRow->chroName = (char *) malloc (sizeof (char) * len + 1);

  strcpy (newRow->chroName, chroName);


  int mls = atoi (locMapLeftStart);
  int mrs = atoi (locMapRightStart);

  if (mls > mrs)		//swap
    {
      newRow->locMapLeftEnd = atoi (locMapRightEnd);
      newRow->locMapRightStart = atoi (locMapLeftStart);
      newRow->locMapLeftStart = atoi (locMapRightStart);
      newRow->locMapRightEnd = atoi (locMapLeftEnd);
      newRow->orientationLeft = orientationRight[0];
      newRow->orientationRight = orientationLeft[0];
    }
  else
    {
      newRow->locMapLeftEnd = atoi (locMapLeftEnd);
      newRow->locMapRightStart = atoi (locMapRightStart);
      newRow->locMapLeftStart = atoi (locMapLeftStart);
      newRow->locMapRightEnd = atoi (locMapRightEnd);
      newRow->orientationLeft = orientationLeft[0];
      newRow->orientationRight = orientationRight[0];
    }

  newRow->avgQual = atof (avgQual);
  newRow->editDistance = atof (editDistance);
  newRow->phredScore = atof (phredScore);

  newRow->libInfo = libInfo;

  len = strlen (readName);
  struct ReadName *r =
    vh_addReadName (hash, readName, newRow->editDistance, newRow->phredScore);
  newRow->readName = r;

  newRow->divetRowId = id;

  return newRow;

}

struct DivetRow *vh_loadDivetRowFromString (struct ReadName *hash[], char *line,
					    struct LibraryInfo *libInfo, int id)
{
  char *readName = strtok (line, DIVET_ROW_DELIMITERS);
  char *chroName = strtok (NULL, DIVET_ROW_DELIMITERS);

  char *locMapLeftStart = strtok (NULL, DIVET_ROW_DELIMITERS);	//skip
  char *locMapLeftEnd = strtok (NULL, DIVET_ROW_DELIMITERS);
  char *orientationLeft = strtok (NULL, DIVET_ROW_DELIMITERS);

  char *chroName2 = strtok (NULL, DIVET_ROW_DELIMITERS);	// This is the same as chroName for IDVE (indicated by =) and for T it is another chr name (We are not deadling with it here)

  char *locMapRightStart = strtok (NULL, DIVET_ROW_DELIMITERS);
  char *locMapRightEnd = strtok (NULL, DIVET_ROW_DELIMITERS);	//skip
  char *orientationRight = strtok (NULL, DIVET_ROW_DELIMITERS);

  char *svType = strtok (NULL, DIVET_ROW_DELIMITERS);	//skip

  char *editDistance = strtok (NULL, DIVET_ROW_DELIMITERS);
  char *avgQual = strtok (NULL, DIVET_ROW_DELIMITERS);
  char *phredScore = strtok (NULL, DIVET_ROW_DELIMITERS);


  struct DivetRow *newRow = createDivetRow (hash, readName, chroName,
					    locMapLeftStart, locMapLeftEnd,
					    orientationLeft,
					    locMapRightStart, locMapRightEnd,
					    orientationRight,
					    svType,
					    editDistance, avgQual, phredScore,
					    libInfo, id);


  return newRow;
}

void vh_freeDivets (struct LibraryInfo *libInfo)
{
  //TODO: To be implemented - FreeDivets
  struct DivetRow *temp = libInfo->head;
  while (temp != NULL)
    {
      //Don't free the readName, as it is shared between more than one divet row
      free (temp->chroName);
      libInfo->head = temp->next;
      free (temp);
      temp = libInfo->head;
    }
  libInfo->head = NULL;
  libInfo->tail = NULL;
  libInfo->size = 0;
}

/**
 * Divet file is read and loaded into the linked lists
 * After running this method, the values of g_headDivet, g_tailDivet and g_divetRowLinkedListSize are set
 */
DivetRow *vh_loadDivetFile (LibraryInfo * libInfo)
{
  FILE *divetFile = safe_fopen (libInfo->libFileAdrs, "r");
  int i;
  
  if (divetFile == NULL)
    {
      sprintf (g_error_message, "Divet file '%s' could not be opened!",
	       libInfo->libFileAdrs);
      vh_quitProgram (EXIT_CODE_DIVET_ERROR);
    }

  char line[400];
  int counter = 0;
  char *token;

  //Initializing the linked list of divet rows
  libInfo->head = NULL;
  libInfo->tail = NULL;
  libInfo->size = 0;

  int counterDivetRow = 0;
  while (!feof (divetFile))
    {

      if (fgets (line, 400, divetFile) > 0)
	;

      if (line == NULL)
	continue;
      // If the read line is empty    
      int isEmpty = 1;
      int len = strlen (line);
      for (i = 0; i < len; i++)
	{
	  if (!isspace (line[i]))
	    {
	      isEmpty = 0;
	      break;
	    }
	}
      if (isEmpty)
	continue;

      DivetRow *newRow = vh_loadDivetRowFromString (libInfo->hash, line, libInfo, counterDivetRow);

      //fixOrientation(newRow);

      if (vh_notInRepeat (newRow) == 1)
	{
	  if (libInfo->head == NULL || libInfo->tail == NULL)
	    {
	      libInfo->head = newRow;
	      libInfo->tail = newRow;
	    }
	  else
	    {
	      libInfo->tail->next = newRow;
	      libInfo->tail = newRow;
	    }
	  libInfo->size++;
	  counterDivetRow++;
	}
      else
	free (newRow);

      line[0] = '\0';

    }

  vh_logInfo ("Finished reading the file");
  fclose (divetFile);
  return libInfo->head;
}

void vh_printDivet (DivetRow * divetRow)
{
  char msg[200];
  sprintf (msg, "%s %s %d %d %c %d %d %c %f %f %g\n",
	   divetRow->readName->readName,
	   divetRow->chroName,
	   divetRow->locMapLeftEnd,
	   divetRow->locMapLeftStart,
	   divetRow->orientationLeft,
	   divetRow->locMapRightStart,
	   divetRow->locMapRightEnd,
	   divetRow->orientationRight,
	   divetRow->editDistance, divetRow->avgQual, divetRow->phredScore);

  vh_logOutput (msg);
}

