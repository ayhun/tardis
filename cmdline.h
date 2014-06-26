#ifndef __COMMANDLINE
#define __COMMANDLINE

int parse_command_line( int, char**, parameters*);
void parse_bam_list( parameters** params);
void tokenize_bam_files( parameters** params);
void print_help( void);

#endif
