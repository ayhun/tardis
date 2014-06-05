#ifndef VH_LOGGER_H__
#define VH_LOGGER_H__

#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#define LOG_LEVEL_NO_LOG 0
#define LOG_LEVEL_LOG_OUTPUT 1
#define LOG_LEVEL_LOG_ERROR 2
#define LOG_LEVEL_LOG_WARNING 4
#define LOG_LEVEL_LOG_INFO 8
#define LOG_LEVEL_DEBUG_INFO 16
#define LOG_LEVEL_ALL 31 

//Defined in vh_logger.c
extern FILE* g_logOutputFile;	//Can be stdout or a log file
//Defined in vh_logger.c
extern int g_currentLogLevel;	//Indicates the current log level that is used in the program
//Defined in vh_logger.c
extern char g_loggerMsgBuffer[];

void logTime();						//Logs the current time
void logOutput(char* message);		//Logs the message without any change
void logError(char* message);		//Logs as: ****ERROR: Message
void logWarning(char* message);		//Logs as: WARN: Message
void logInfo(char* message);		//Logs as: INFO: Message
void logDebug(char* message);		//Logs as: DEBUG: Message
void logDebugint(int i);		//Logs as: DEBUG: Message
void log(char* message, int logLevel);
void initLogger(FILE* logOutputFile, int logLevel);

#endif

