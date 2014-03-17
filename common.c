#include <stdio.h>
#include <stdlib.h>
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
}

void print_params(parameters *params){
  printf("bamFile: %s\n", params->bamFile);
  printf("refGenome: %s\n", params->refGenome);
  printf("reps: %s\n", params->reps);
  printf("dups: %s\n", params->dups);
  printf("gaps: %s\n", params->gaps);
  printf("mei: %s\n", params->mei);
}
