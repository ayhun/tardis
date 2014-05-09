#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"

void load_config()
{
	FILE* config;
	char* next_line = NULL;
	char config_filename[MAX_LENGTH];
	int len = 0;
	int bytes_read;
	int i;
	
	/* Initialize the executable paths to empty */
	path_samtools[0] = 0;
	path_bcftools[0] = 0;
	path_mrfast[0] = 0;

	/* Combine the home directory path with the name of the configuration file */
	sprintf( config_filename, "%s/%s", getenv( "HOME"), CONFIG_FILE);

	/* Open the configuration file for reading */
	config = fopen( config_filename, "r");
	if( config == NULL)
	{
		/* Create new config file */
		create_config( config_filename);
	}
	else
	{
		/* Get the paths from the pre-existing configuration file */
		i = 0;
		while( ( bytes_read = getline( &next_line, &len, config)) != -1)
		{
			if( next_line[0] != '#')
			{
				if( i == 0)
				{
					strcpy( path_samtools, next_line);
				}
				else if( i == 1)
				{
					strcpy( path_bcftools, next_line);
				}
				else if( i == 2)
				{
					strcpy( path_mrfast, next_line);
				}
				else
				{
					fprintf( stderr, "Configuration file has wrong format or unspecified external tools used.\n");
				}
				i = i + 1;
			} 
		}

		/* Free memory allocated internally by getline */
		free( next_line);

		/* If the first characters of the paths are still '0', then they are either not
		   installed, or not in the PATH. */
		if( path_samtools[0] == 0)
		{
			fprintf( stderr, "Warning: samtools path is not in the configuration file.\n");
		}
		else
		{
			fprintf( stderr, "samtools path: %s\n", path_samtools);
		}

		if( path_bcftools[0] == 0)
		{
			fprintf( stderr, "Warning: bcftools path is not in the configuration file.\n");
		}
		else
		{
			fprintf( stderr, "bcftools path: %s\n", path_bcftools);
		}

		if( path_mrfast[0] == 0)
		{
			fprintf( stderr, "Warning: mrfast path is not in the configuration file.\n");
		}
		else
		{
			fprintf( stderr, "mrfast path: %s\n", path_mrfast);
		}
	}
}

void create_config( char* config_filename)
{
	FILE* config;
	FILE* pipe;
  	char executable_path[MAX_LENGTH];

	/* popen(...) executes shell commands */
	/* "which" finds the path to the specified executable if it exists */
	/* "2>/dev/null" redirects the output into a virtual black hole */
	pipe = popen( "which samtools 2>/dev/null", "r");
	if( pipe == NULL)
	{
		fprintf( stderr, "Error opening pipe\n");
	}
	else
	{
		if( fgets( executable_path, MAX_LENGTH, pipe) == NULL)
		{
			fprintf( stderr, "samtools not found in PATH. Install it or manually configure the %s file.\n", config_filename);
		}
		else
		{
			strcpy( path_samtools, executable_path);
		}
		pclose( pipe);
	}

	pipe = popen( "which bcftools 2>/dev/null", "r");
	if( pipe == NULL)
	{
		fprintf( stderr, "Error opening pipe\n");
	}
	else
	{
		if( fgets( executable_path, MAX_LENGTH, pipe) == NULL)
		{
			fprintf( stderr, "bcftools not found in PATH. Install it or manually configure the %s file.\n", config_filename);
		}
		else
		{
			strcpy( path_bcftools, executable_path);
		}
		pclose( pipe);
	}

	pipe = popen( "which mrfast 2>/dev/null", "r");
	if( pipe == NULL)
	{
		fprintf( stderr, "Error opening pipe\n");
	}
	else
	{
		if( fgets( executable_path, MAX_LENGTH, pipe) == NULL)
		{
			fprintf( stderr, "mrfast not found in PATH. Install it or manually configure the %s file.\n", config_filename);
		}
		else
		{
			strcpy( path_mrfast, executable_path);
		}
		pclose( pipe);
	}
	
	config = fopen( config_filename, "w");
	fprintf( config, "SAMTOOLS = %s", path_samtools);
	fprintf( config, "BCFTOOLS = %s", path_bcftools);
	fprintf( config, "MRFAST = %s", path_mrfast);
	fclose( config);

	fprintf( stderr,"\n");
	fprintf( stderr, "*************************************************************************\n");
	fprintf( stderr, "*           Config file $HOME/.tardis_config is created.                *\n"); 
	fprintf( stderr, "*            Check the configuration file for any errors,               *\n");
	fprintf( stderr, "* then run TARDIS again to load the configuration and normal operation. *\n");
	fprintf( stderr, "*************************************************************************\n");
}
