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
	char svfile[2000];
	for( i = 0; i < params->num_bams; i++)
	{
		
		FILE *divetlib;
		divetlib=fopen(divetadd,"w");
		fprintf( divetlib, "%d \n", in_bams[i]->num_libraries);
		for ( j = 0; j < in_bams[i]->num_libraries; j++){
			fprintf(divetlib, "%s %s %s-%s.sam_DIVET.vh %d %d %d\n",in_bams[i]->libraries[j]->libname,in_bams[i]->sample_name,in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname,in_bams[i]->libraries[j]->conc_min, in_bams[i]->libraries[j]->conc_max, in_bams[i]->libraries[j]->read_length);
			}
			fclose(divetlib);
			for ( j = 0; j < in_bams[i]->num_libraries; j++){
			  sprintf(outputfile,"%s-%s.cluster.out",in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname);
			  sprintf(outputread,"%s-%s.name",in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname);
			  sprintf(svfile,"%s-%s.out.sv",in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname);
			  vh_clustering (divetadd, in_bams[i], params->gaps, "build37.satellite.bed", preProsPrune, outputfile, outputread, overMapLimit);
			  vh_setcover(divetadd, outputread, outputfile, svfile);
			}
			remove(divetadd);
	}  
	
	return RETURN_SUCCESS;
}
