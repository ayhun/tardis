#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vh.h"

int run_vh(parameters *params, bam_info ** in_bams){
	int i, j;
	double preProsPrune= 0.001;
	int overMapLimit=500;
	int return_value;
	char divetfile[MAX_SEQ];
	char outputfile[MAX_SEQ];
	char outputread[MAX_SEQ];
	char svfile[MAX_SEQ];

	for( i = 0; i < params->num_bams; i++)
	{
		
		for ( j = 0; j < in_bams[i]->num_libraries; j++)
		{
		        sprintf(divetfile, "%s-%s.sam_DIVET.vh",  in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname);
			set_str(&(in_bams[i]->libraries[j]->divet), divetfile);
		}
		for ( j = 0; j < in_bams[i]->num_libraries; j++)
		{
		  fprintf(stderr, "Calculating maximal clusters for sample: %s, library: %s.\n", in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname);
		  sprintf(outputfile,"%s-%s.cluster.out",in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname);
		  sprintf(outputread,"%s-%s.name",in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname);
		  sprintf(svfile,"%s-%s.out.sv",in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname);
		  if ( !params->skip_vhcluster)
		    vh_clustering (in_bams[i], params->gaps, params->reps, preProsPrune, outputfile, outputread, overMapLimit);
		}
		//remove(outputfile);
		//remove(outputread);
	}  
	
	fprintf(stderr, "Applying SET-COVER to find putative structural variation.\n");
	vh_setcover(in_bams, params->num_bams, outputread, outputfile, svfile); 
	fprintf(stderr, "VariationHunter run is complete.\n");
	
	return RETURN_SUCCESS;
}
