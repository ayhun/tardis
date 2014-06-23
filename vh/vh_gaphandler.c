#include "stdio.h"
#include "vh_common.h"
#include "vh_gaphandler.h"


Gap g_gapTable[1000];		//TODO: Malloc instead
int g_gapTableSize = 0;

Chro g_chroTable[100];
int g_chroTableSize = 0;




void vh_readGapTable (char *gapFileName)
{


  FILE *gapFile = fopen (gapFileName, "r");
  if (!gapFile)
    //TODO: handle me
    ;

  int index = 0;
  while (!feof (gapFile))
    {
      //int res = fscanf(gapFile,"%s%d%d", chroName, &(g_gapTable[index].start), &(g_gapTable[index].end));

      int res = fscanf (gapFile, "%s%d%d\n", g_gapTable[index].chroName,
			&(g_gapTable[index].start), &(g_gapTable[index].end));
      //g_gapTable[index].chroName = (char *) malloc((strlen(chroName)+1)*sizeof(char));
      //strcpy(g_gapTable[index].chroName, chroName);
//              printf("GapTable: %s %i %i\n",g_gapTable[index].chroName, g_gapTable[index].start, g_gapTable[index].end);
      if (res != 3)
	break;
      //logDebug(g_gapTable[index].chroName); 
      //g_chroTable[index].chroName = (char *) malloc((strlen(chrName)+1)*sizeof(char));
      //strcpy(g_chroTable[index].chroName, chrName);

      index++;
    }
  g_gapTableSize = index;

  fclose (gapFile);
}




void vh_readInitFile (char *fileName)
{
  FILE *initFile = fopen (fileName, "r");

  while (!feof (initFile))
    {
      int i = 0;
      char str[50];
      int res = fscanf (initFile, "%s%d", str, &i);
      if (res != 2)
	break;

      if (!strcmp (str, "max_potential_cluster_size"))	//TODO: Move as a constant to gap_handler.c
	g_maxListBrkPointIntr = i;
    }

  fclose (initFile);
}

void vh_readChros (char *fileName)
{
//      char chroName[500];

  FILE *chroFile = fopen (fileName, "r");
  if (!chroFile)
    //TODO: handle me
    ;

  int index = 0;
  while (!feof (chroFile))
    {
      int res = fscanf (chroFile, "%s%d", g_chroTable[index].chroName,
			&(g_chroTable[index].size));
      //int res = fscanf(chroFile,"%s%d", chroName, &(g_chroTable[index].size));
      //g_chroTable[index].chroName = (char *) malloc((strlen(chroName)+1)*sizeof(char));
      //strcpy(g_chroTable[index].chroName, chroName);
      if (res != 2)		//TODO: Error in file
	break;
      //logDebug(g_chroTable[index].chroName);
      index++;
    }
  g_chroTableSize = index;

  fclose (chroFile);
}
