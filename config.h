#ifndef __CONFIG
#define __CONFIG

/* Name of the configuration file can be changed here */
#define CONFIG_FILE ".tardis_config" 

/* Maximum length of a single line of the configuration file */
#define MAX_LENGTH 1024

/* External tool executable paths */
char path_samtools[MAX_LENGTH];
char path_bcftools[MAX_LENGTH];
char path_mrfast[MAX_LENGTH];

/* Function Prototypes */
void load_config();
void create_config( char* config_filename);

#endif
