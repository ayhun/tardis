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

	sprintf(outputfile,"%s.clusters", params->outprefix);
	sprintf(outputread,"%s.name", params->outprefix);
	sprintf(svfile,"%s.sv", params->outprefix);

	for( i = 0; i < params->num_bams; i++)
	{
		
		for ( j = 0; j < in_bams[i]->num_libraries; j++)
		{
		        sprintf(divetfile, "%s-%s.sam_DIVET.vh",  in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname);
			set_str(&(in_bams[i]->libraries[j]->divet), divetfile);
		}
	}  

	fprintf( stderr, "Now running VariationHunter/CommonLAW...\n");	
	fprintf( stderr, "Calculating maximal clusters.\n");
	if ( TARDIS_DEBUG && !params->skip_vhcluster) // this parameter is only intended for debugging purposes. End users shouldn't use this
	  vh_clustering (in_bams, params->num_bams, params->gaps, params->reps, preProsPrune, outputfile, outputread, overMapLimit);
	fprintf(stderr, "Applying SET-COVER approximation to find putative structural variation.\n");
	vh_setcover(in_bams, params->num_bams, outputread, outputfile, svfile); 

	fprintf(stderr, "VariationHunter/CommonLAW run is complete. Results are in the %s file.\n", svfile);

	if (! TARDIS_DEBUG)
	{
	  remove(outputfile);
	  remove(outputread);
	}
	
	return RETURN_SUCCESS;
}
