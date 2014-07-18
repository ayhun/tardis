#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "clustering.h"
int vhclustering(parameters *params, bam_info ** in_bams){
	int i, j;
	double preProsPrune= 0.001;
	int overMapLimit=500;
	int return_value;
	char* divetadd= "divetlib";
	char outputfile[2000];
	char outputread[2000];
	for( i = 0; i < params->num_bams; i++)
	{
		
		int read_length = 101;
		FILE *divetlib;
		divetlib=fopen(divetadd,"w");
		fprintf( divetlib, "%d \n", in_bams[i]->num_libraries);
		for ( j = 0; j < in_bams[i]->num_libraries; j++){
			fprintf(divetlib, "%s %s %s-%s.sam_divet.vh %d %d %d\n",in_bams[i]->libraries[j]->libname,in_bams[i]->sample_name,in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname,in_bams[i]->libraries[j]->conc_min, in_bams[i]->libraries[j]->conc_max, read_length);
			}
			fclose(divetlib);
			for ( j = 0; j < in_bams[i]->num_libraries; j++){
			sprintf(outputfile,"%s-%s.cluster.out",in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname);
			sprintf(outputread,"%s-%s.name",in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname);
			  printf("check 1\n");
			vh_clustering (divetadd, in_bams[i], params->gaps, params->reps, preProsPrune, outputfile, outputread, overMapLimit);
			vh_setcover(divetadd, outputread, outputfile);
		}
			remove(divetadd);
	}  
	
	return RETURN_SUCCESS;
}