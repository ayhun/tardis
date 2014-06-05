#include "vh_logger.h"

char tempMsg[500];
FILE* g_logOutputFile =stdout;
//C++ to C conversion error
// g_logOutputFile=stdout;
int g_currentLogLevel = LOG_LEVEL_ALL;
char g_loggerMsgBuffer[400];

void log(char* message, int logLevel)
{
	if (!(logLevel & g_currentLogLevel))
		return;

	fprintf(g_logOutputFile, message);
}

void initLogger(FILE* logOutputFile, int logLevel)
{
	g_logOutputFile = logOutputFile;
	g_currentLogLevel = logLevel;
}

void logError(char* message)
{
	if (!(LOG_LEVEL_LOG_ERROR & g_currentLogLevel))
		return;

	sprintf(tempMsg, "****ERROR: \t%s\n", message);
	log(tempMsg, LOG_LEVEL_LOG_ERROR);
}

void logOutput(char* message)
{
	log(message, LOG_LEVEL_LOG_OUTPUT);
}

void logInfo(char* message)
{
	if (!(LOG_LEVEL_LOG_INFO & g_currentLogLevel))
		return;

	sprintf(tempMsg, "INFO: \t%s\n", message);
	log(tempMsg, LOG_LEVEL_LOG_INFO);
}

void logWarning(char* message)
{
	if (!(LOG_LEVEL_LOG_WARNING & g_currentLogLevel))
		return;

	sprintf(tempMsg, "WARN: \t%s\n", message);
	log(tempMsg, LOG_LEVEL_LOG_WARNING);
}

void logDebugint(int i)
{
	if (!(LOG_LEVEL_DEBUG_INFO& g_currentLogLevel))
		return;

	sprintf(tempMsg, "DEBUG \t%d\n", i);
	log(tempMsg, LOG_LEVEL_DEBUG_INFO);
}

void logDebug(char* message)
{
	if (!(LOG_LEVEL_DEBUG_INFO& g_currentLogLevel))
		return;

	sprintf(tempMsg, "DEBUG \t%s\n", message);
	log(tempMsg, LOG_LEVEL_DEBUG_INFO);
}

void logTime()
{
	log("Log time to be implemented", LOG_LEVEL_LOG_WARNING);
}

