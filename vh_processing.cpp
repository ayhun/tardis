#include <string.h>
#include "vh_divethandler.h"
#include "vh_maximalCluster.h"
#include "vh_gaphandler.h"
#include "vh_repeathandler.h"
#include "vh_processing.h"

char **allReadNameList;

void printHelp()
{
	FILE* helpFile = fopen("vh_help.hlp", "r");
	if (helpFile == NULL)
	{
		logError("Help file \"vh_help.hlp\" not found.\n");
	}
	else
	{
		char line[1024];
		while (!feof(helpFile))
		{
			line[0] = line[1] = '\0';
			fgets(line, 1024, helpFile);
			logOutput(line);
		}
		fclose(helpFile);
	}
}
LibraryInfo* readLibraryInfos(char* libFileAdrs)
{
	FILE* libFile = fopen(libFileAdrs, "r");
	if (libFile == NULL)
	{
		sprintf(g_error_message, "Lib file cannot be opened: %s", libFileAdrs);
		quitProgram(EXIT_CODE_DIVET_ERROR);
		return NULL;
	}

	int NUM_OF_FIELDS = 6;
	int numOfLines = 0;
	fscanf(libFile, "%d", &numOfLines);
	for (int lineNumber = 0; lineNumber < numOfLines; lineNumber++) {
		char nameStr[100];
		char indName[100];
		char adrsStr[200];
		int min, max, readLen;
		int numRead = fscanf(libFile, "%s %s %s %d %d %d", nameStr, indName, adrsStr, &min, &max, &readLen);
		if (numRead == NUM_OF_FIELDS)
		{
			FILE* testfile = fopen(adrsStr, "r");
			if (testfile == NULL)
			{
				sprintf(g_error_message, "File %s cannot be opened!", adrsStr);
				quitProgram(EXIT_CODE_DIVET_ERROR);
			} else {
				fclose(testfile);
			}
			//create new LibInfo object
			LibraryInfo* newLibInfo = (LibraryInfo*)malloc(sizeof(LibraryInfo));
			strcpy(newLibInfo->libName, nameStr);
			strcpy(newLibInfo->libFileAdrs, adrsStr);
			printf("%s %s\n", newLibInfo->libName, newLibInfo->libFileAdrs);
			newLibInfo->minDelta = min-2*readLen;
			newLibInfo->maxDelta = max-2*readLen;
			if (newLibInfo->minDelta<0)
				newLibInfo->minDelta=0;
			if (newLibInfo->maxDelta<0)
				newLibInfo->maxDelta=0;
			newLibInfo->readLen = readLen;
			newLibInfo->hash = (ReadName**)malloc(NHASH * sizeof(ReadName*));
			for (int i = 0; i < NHASH; i++)
				newLibInfo->hash[i] = NULL;
			newLibInfo->head = NULL;
			newLibInfo->tail = NULL;
			newLibInfo->size = 0;
			newLibInfo->next = NULL;

			//If first one in the linkedlist
			if (g_libInfo == NULL)
			{
				g_libInfo = newLibInfo;
			}
			else  //add to the end of the linked list
			{
				LibraryInfo* t;
				for (t = g_libInfo; t->next != NULL; t=t->next); //Skip till end of LinkedList

				t->next = newLibInfo;
			}
		}
		if (feof(libFile))
			break;
		if (numRead != NUM_OF_FIELDS) //Error in reading from file
		{
			sscanf(g_error_message, "Error in lib file - around line: %d", lineNumber);
			logError(g_error_message);
			quitProgram(EXIT_CODE_DIVET_ERROR);
			break;
		}
	}

	return g_libInfo;
}
void printVersion()
{
	logOutput("Variation Hunter 3.0 New Version\r\n");
}

void quitProgram(int exitCode)
{
	if (exitCode != 0)
	{
		if (!strcmp(g_error_message, ""))
			strcpy(g_error_message, "Unexpected Error Found!\n");

		char errMsg[150];
		sprintf(errMsg, "Program Exited with Errors:%d", exitCode); 
		logError(errMsg);
		switch (exitCode)
		{
			case EXIT_CODE_ARG_ERROR:
				sprintf(errMsg, "Error in options: %s", g_error_message);
				logError(errMsg);
				printHelp();
				break;
			case EXIT_CODE_DIVET_ERROR:
				sprintf(errMsg, "Error in DIVET file: %s", g_error_message);
				logError(errMsg);
				printHelp();
				break;
			case EXIT_CODE_MEMORY_ERROR:
				sprintf(errMsg, "Memory Problem Occured: %s", g_error_message);
				logError(errMsg);
				break;
			default:
				sprintf(errMsg, "Uncategorized Error Found: %s", g_error_message);
				logError(errMsg);
				printHelp();
				break;
		}
		exit(exitCode);
	}
	else 
	{
		logInfo("Program Exited Successfully.\r\n");
		exit(0);
	}
}

//Gets the adrs to a library info file,reads the file and loads the contents into a linkedlist of LibraryInfo objects
//ALSO, sets the g_libInfo reference 


void pruneAndNormalizeDivets(LibraryInfo* lib, double preProsPrune, int overMapLimit)
{
	DivetRow* cursor = lib->head->next;

	

	if (cursor == NULL)
		return;
	printf("%f\n", preProsPrune);


	//We are not checking the first and last one
	while (cursor != NULL && cursor->next != NULL)
	{
		if ( ((cursor->next->phredScore / cursor->next->readName->sumPhredValue) < preProsPrune)||(cursor->next->readName->occurrences > overMapLimit))
		{
			DivetRow* toBeDeleted = cursor->next;
			cursor->next = cursor->next->next;
			lib->size--;	
			free(toBeDeleted->chroName);
			free(toBeDeleted);
		} 
		else 
		{
			cursor->next->phredScore = cursor->next->phredScore / cursor->next->readName->sumPhredValue;	
			cursor = cursor->next;
		}
			
		
	}
}


static int cmprReadNameStr(const void *a, const void *b) {
	return strcmp(*(char **)a, *(char **)b);
}



//We assume that the arguments are valid at this point and the file exists
void run(char* libFileAdrs, char* chroFileName, char* gapFileName, char * repeatFileName,char* initFileName, double preProsPrune, double minSVSup, char* outputFile, char *outputRead, int overMapLimit)
{
	int totalNumUniqueReads=0;
	int indexStart=0;
	readInitFile(initFileName);
	readChros(chroFileName);
	readGapTable(gapFileName);
	readRepeatTable(repeatFileName);
	readLibraryInfos(libFileAdrs);
	LibraryInfo* cursor = g_libInfo;
	fileOutput=fopen(outputFile,"w");
	printf("prune %f\n", preProsPrune); 
	for (;cursor;cursor = cursor->next)
	{
		logInfo("Reading Divet Fileshyus ...");
		loadDivetFile(cursor);
		sprintf(g_loggerMsgBuffer, "%d rows loaded successfully.", cursor->size);
		logInfo(g_loggerMsgBuffer);
		pruneAndNormalizeDivets(cursor, preProsPrune, overMapLimit);
		sprintf(g_loggerMsgBuffer, "%d rows after pruning.", cursor->size);
		logInfo(g_loggerMsgBuffer);
	}
	
	if (strcmp(outputRead,""))
	{
		logInfo("Writing ReadName Sorted");
		cursor=g_libInfo;
		FILE *fileOutputReadName;
		fileOutputReadName=fopen(outputRead,"w");	
		while(cursor!=NULL)
		{
			totalNumUniqueReads=totalNumUniqueReads+countNumReads(cursor->hash);
			cursor=cursor->next;
		}
		cursor=g_libInfo;
		allReadNameList=(char **) malloc(totalNumUniqueReads * sizeof(char *));
		for(;cursor;cursor = cursor->next)
		{
			indexStart=exportToArray(cursor->hash, allReadNameList, indexStart);
		}
		qsort(allReadNameList, totalNumUniqueReads, sizeof(char *), cmprReadNameStr);
		fprintf(fileOutputReadName, "%i\n", totalNumUniqueReads+1);
		for (int count=0; count<totalNumUniqueReads; count++)
		{
			fprintf(fileOutputReadName, "%s\n", allReadNameList[count]);
		}
		fclose(fileOutputReadName);
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
	for (int i = 0; i < g_chroTableSize; i++)
	{
		initializeReadMapping_Deletion(g_chroTable[i].chroName, g_chroTable[i].size);
		createDeletionClusters(g_chroTable[i].size);
		finalizeReadMapping(g_chroTable[i].chroName, g_chroTable[i].size);
		initializeReadMapping_Inversion(g_chroTable[i].chroName, g_chroTable[i].size);
		createInversionClusters(g_chroTable[i].size);
		finalizeReadMapping(g_chroTable[i].chroName, g_chroTable[i].size);
		initializeReadMapping_Insertion(g_chroTable[i].chroName, g_chroTable[i].size);
		createInsertionClusters(g_chroTable[i].size);
		finalizeReadMapping(g_chroTable[i].chroName, g_chroTable[i].size);
	}
	fclose(fileOutput);

}

bool isValid(MainOptions o)
{
	//TODO: To be completed. Check for files to exists, etc.
	if (!strcmp(o.libFileAdrs, "") || !strcmp(o.initializeFileName, "") || !strcmp(o.chroFileName, ""))
		return false;
	return true;
}

#ifdef MAIN_MAIN
int vhprocessing(int argc, char** argv)
{
	initLogger(stdout, LOG_LEVEL_ALL);

	//overMapLimit = 1000; // default value. The user can change this value.

	MainOptions mainOptions = {0, 0, "", "", "", "", "", "", "", false, false};

	int errorCode  ;
	// = parse_command_line( argc, argv, &params, &mainOptions);
{
	if (errorCode != 0)
		quitProgram(EXIT_CODE_ARG_ERROR); 
	else if (mainOptions.helpWanted)
		printHelp();
	else if (mainOptions.versionWanted)
		printVersion();
	else 
	{
		if (!isValid(mainOptions))
		{
			strcpy(g_error_message, "Error in arguments!");
			quitProgram(EXIT_CODE_ARG_ERROR);
		}
		//TODO: Test the arguments are OK and the input is valid
		run(mainOptions.libFileAdrs, mainOptions.chroFileName, mainOptions.gapFileName, mainOptions.repeatFileName, mainOptions.initializeFileName, mainOptions.prunProb, mainOptions.svSup, mainOptions.outputFile, mainOptions.outputRead, mainOptions.overMapLimit);
		quitProgram(EXIT_CODE_SUCCESS);
	}
}
}
#endif
