
#ifndef VH_DIVET_HANDLER__
#define VH_DIVET_HANDLER__

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>		//For isspace

#include "vh_common.h"
#include "vh_hash.h"

#define DIVET_ROW_DELIMITERS " \t\r\n"
#define BUFFER_SIZE 20000

#define MAX_LIB_NAME_LEN  200
#define MAX_FILE_NAME_LEN 200

typedef struct LibraryInfo
{
  char libName[MAX_LIB_NAME_LEN];
  char libFileAdrs[MAX_FILE_NAME_LEN];
  int minDelta;
  int maxDelta;
  int readLen;
  int size;

  struct ReadName **hash;
  struct DivetRow *head;
  struct DivetRow *tail;

  struct LibraryInfo *next;
} LibraryInfo;


typedef struct DivetRow
{
  struct ReadName *readName;

  char *chroName;
  int locMapLeftEnd;
  int locMapLeftStart;
  char orientationLeft;
  int locMapRightStart;
  int locMapRightEnd;
  char orientationRight;

  double avgQual;
  double editDistance;
  double phredScore;
  ////////////////////////I added this variable to give an
  int divetRowId;

  struct LibraryInfo *libInfo;

  struct DivetRow *next;
} DivetRow;


struct DivetRow *createDivetRow (struct ReadName *hash[],
				 char *readName,
				 char *chroName,
				 char *locMapLeftStart,
				 char *locMapLeftEnd,
				 char *orientatinoLeft,
				 char *locMapRightStart,
				 char *locMapRightEnd,
				 char *orientationRight,
				 char *svType,
				 char *editDistance,
				 char *avgQual,
				 char *phredScore,
				 struct LibraryInfo *libInfo, int id);

struct DivetRow *vh_loadDivetRowFromString (struct ReadName *hash[], char *line,
					 struct LibraryInfo *libInfo, int id);
void vh_freeDivets ();
struct DivetRow *vh_loadDivetFile (struct LibraryInfo *libInfo);
void vh_printDivet (struct DivetRow *divetRow);
//void fixOrientation(DivetRow* row);

#endif
