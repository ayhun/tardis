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
  if (!gapFile);
    //TODO: handle me
  

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




void vh_readInitFile ()
{
  // FILE *initFile = fopen (fileName, "r");
  int i = 10000000;
 //  while (!feof (initFile))
 //    {
 //      int i = 0;
 //      char str[50];
 //      int res = fscanf (initFile, "%s%d", str, &i);
 //      if (res != 2)
	// break;

 //      if (!strcmp (str, "max_potential_cluster_size"))	//TODO: Move as a constant to gap_handler.c
	g_maxListBrkPointIntr = i;
    // }

  // fclose (initFile);
}

void vh_readChros (bam_info* in_bam){
  int i;
  for(i=0;i<in_bam->num_chrom;i++){
    strcpy(g_chroTable[i].chroName,in_bam->chrom_names[i]);
    g_chroTable[i].size=in_bam->chrom_lengths[i];
  }
  g_chroTableSize=i;

}
