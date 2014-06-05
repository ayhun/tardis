
#ifndef VH_DIVET_HANDLER__
#define VH_DIVET_HANDLER__

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h> //For isspace
#include <stdbool.h>
#include "vh_common.h"
#include "vh_hash.h"

#define DIVET_ROW_DELIMITERS " \t\r\n"
#define BUFFER_SIZE 20000

#define MAX_LIB_NAME_LEN  200
#define MAX_FILE_NAME_LEN 200

struct LibraryInfo {
	char libName[MAX_LIB_NAME_LEN];
	char libFileAdrs[MAX_FILE_NAME_LEN];
	int minDelta;
	int maxDelta;
	int readLen;

	ReadName** hash;
	DivetRow* head;
	DivetRow* tail;
	int size;

	LibraryInfo* next;
}; 

struct LibraryInfo;

struct DivetRow
{
	ReadName* readName;

	char* chroName;
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

	LibraryInfo* libInfo;

	struct DivetRow* next;
};


DivetRow* createDivetRow(
	ReadName* hash[],
	char* readName,
	char* chroName,
	char* 	locMapLeftStart,
	char* 	locMapLeftEnd,
	char* 	orientatinoLeft,
	char* 	locMapRightStart,
	char* 	locMapRightEnd,
	char* 	orientationRight,
	char* 	svType,
	char* 	editDistance,	
	char* 	avgQual, 	
	char* 	phredScore,
	int id);

DivetRow* loadDivetRowFromString(ReadName* hash[], char* line, int id);
void freeDivets();
DivetRow* loadDivetFile(LibraryInfo* libInfo);
void printDivet(DivetRow* divetRow);
//void fixOrientation(DivetRow* row);

#endif

