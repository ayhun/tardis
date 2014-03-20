#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "common.h"

void init_params(parameters **params){
  /* initialize parameters */
  *params = (parameters *) malloc(sizeof(parameters));
  (*params)->refGenome = NULL;
  (*params)->reps = NULL;
  (*params)->dups = NULL;
  (*params)->bamFile = NULL;
  (*params)->gaps = NULL;
  (*params)->mei = NULL;
  (*params)->sampleGender = MALE;
  (*params)->runVH = 0; 
  (*params)->runNS = 0;
  (*params)->runSR = 0;
  (*params)->threads = 1;
}

void print_params(parameters *params){
  printf("bamFile: %s\n", params->bamFile);
  printf("refGenome: %s\n", params->refGenome);
  printf("reps: %s\n", params->reps);
  printf("dups: %s\n", params->dups);
  printf("gaps: %s\n", params->gaps);
  printf("mei: %s\n", params->mei);
}

//void print_error(char *msg){
void print_error(char *msg){
  /*
    print error message and exit
   */

  fprintf(stderr, "\n%s\n", msg);
  fprintf(stderr, "Invoke parameter -h for help.\n");
  exit (0);
}


FILE *gfOpen(char *fname, char *mode){
  /*
    gfOpen: graceful file open.
    try to open a file; exit if file does not exist
   */
  FILE *file;
  char err[500];
  file = fopen(fname, mode);
  
  if (file == NULL){
    sprintf(err, "[TARDIS INPUT ERROR] Unable to open file %s in %s mode.", fname, mode[0]=='w' ? "write" : "read");
    print_error(err);
  }
  return file;
}
