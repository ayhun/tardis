#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include "tardis.h"
#include "vh2vcf.h"

#define MAX 300000000

char repeats[MAX];
char dups[MAX];


float calc_dup(int, int);
float calc_rep(int, int);


char header[20000] = "##fileformat=VCFv4.1\n##ALT=<ID=DEL,Description=\"Deletion\">\n##FILTER=<ID=LowQual,Description=\"Genotype call confidence below LOD 1.3\">\n##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">\n##FORMAT=<ID=CNQ,Number=1,Type=Float,Description=\"Copy number genotype quality for imprecise events\">\n##FORMAT=<ID=FT,Number=.,Type=String,Description=\"Per-sample genotype filter, PASS for called genotypes or list of excluding filters\">\n##FORMAT=<ID=GL,Number=3,Type=Float,Description=\"Genotype Likelihoods\">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##INFO=<ID=BKPTID,Number=.,Type=String,Description=\"ID of the assembled alternate allele in the assembly file\">\n##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">\n##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n##INFO=<ID=DBRIPID,Number=1,Type=String,Description=\"ID of this element in DBRIP\">\n##INFO=<ID=DBVARID,Number=1,Type=String,Description=\"ID of this element in DBVAR\">\n##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample ID\">\n##INFO=<ID=SVALG,Number=1,Type=String,Description=\"SV discovery algorithm\">\n##INFO=<ID=SUP,Number=1,Type=Integer,Description=\"Number of supporting read pairs\">\n##INFO=<ID=DGVID,Number=1,Type=String,Description=\"ID of this element in Database of Genomic Variation\">\n##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate of this variant\">\n##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">\n##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical micro-homology at event breakpoints\">\n##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n##INFO=<ID=MEINFO,Number=4,Type=String,Description=\"Mobile element info of the form NAME\">\n##INFO=<ID=METRANS,Number=4,Type=String,Description=\"Mobile element transduction info of the form CHR\">\n##INFO=<ID=NOVEL,Number=0,Type=Flag,Description=\"Indicates a novel structural variation\">\n##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n##INFO=<ID=SVMETHOD,Number=.,Type=String,Description=\"Type of approach used to detect SV: RP (read pair), RD (read depth), SR (split read), or AS (assembly)\">\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n##INFO=<ID=VALIDATED,Number=0,Type=Flag,Description=\"Validated by PCR, Assembly or Microarray\">\n##INFO=<ID=VALMETHOD,Number=.,Type=String,Description=\"Type of validation: CGH, PCR, SAV (superarray), CAP (capture-array), or ASM (assembly)\">\n##reference=1000GenomesPilot-NCBI37\n";

char header2[1000]= "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";


FILE *divet;
FILE *vcf;
FILE *rep;
FILE *dupfile;

char samplename[100];

int id = 1;
int total_del_bp=0;
int total_ins_bp=0;
int min_del=100000;
int max_del=0;
int del_cnt = 0;
int ins_cnt = 0;
int inv_cnt = 0;
int ram_for_del = 1;
int i,j;



int conversion(parameters *params, bam_info ** in_bams, int i, int j){
	int minsup=2;  //??
	int maxsup=100;  //??
	char fname[200];
	char outfname[200];
	char rname[200];
	char dname[200];
	float maxdup = 0.3;
	int s, e;
	char chrom[50];
	char thischrom[50];
	int os, oe, is, ie;
	char type;
	int support;
	char word[100];
	float dup_content;
	float rep_content_1;
	float rep_content_2;
	float prune = 10000;  //??

	time_t rawtime;
	struct tm * timeinfo;

	int mininvsup=-1;  //??

	time ( &rawtime );
	timeinfo = localtime ( &rawtime );


	strcpy(samplename, "VHcall_");

	sprintf(fname,"%s-%s.sam_DIVET.vh",in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname);
	sprintf(outfname,"%s-%s.vcf",in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname);

	// for (i=1;i<argc;i++){
	//   if (!strcmp(argv[i], "-i")) strcpy(fname, argv[i+1]);
	//   else if (!strcmp(argv[i], "-o")) strcpy(outfname, argv[i+1]);
	//   else if (!strcmp(argv[i], "-r")) strcpy(rname, argv[i+1]);
	//   else if (!strcmp(argv[i], "-d")) strcpy(dname, argv[i+1]);
	//   else if (!strcmp(argv[i], "-min")) minsup = atoi(argv[i+1]);              ??
	//   else if (!strcmp(argv[i], "-minv")) mininvsup = atoi(argv[i+1]);          ??
	//   else if (!strcmp(argv[i], "-max")) maxsup = atoi(argv[i+1]);              ??
	//   else if (!strcmp(argv[i], "-f")) maxdup = atof(argv[i+1]);                ??
	//   else if (!strcmp(argv[i], "-p")) prune = atof(argv[i+1]);				   ??
	//   else if (!strcmp(argv[i], "-xx")) male = 0;
	//   else if (!strcmp(argv[i], "-s")) strcpy(samplename, argv[i+1]);		   ??
	//   else if (!strcmp(argv[i], "-h")) {printHelp(argv[0]); return 0;}
	//   else if (!strcmp(argv[i], "-noram")) ram_for_del = 0;
	// }

    if (mininvsup == -1)
      mininvsup = minsup;
  

    if (fname[0]==0){
      fprintf(stderr, "VH_calls file needed through -i\n"); return -1;
    }
  
      if (outfname[0]==0){
        fprintf(stderr, "Output file needed through -o\n"); return -1;
      }
  
      if (rname[0]==0){
        fprintf(stderr, "Transposons BED file needed through -r\n"); return -1;
      }
  
      if (dname[0]==0){
        fprintf(stderr, "SegDups BED file needed through -d\n"); return -1;
      }
  
  
      divet = fopen (fname, "r");
      vcf = fopen (outfname, "w");
      rep = fopen (params->reps, "r");
      dupfile = fopen (params->dups, "r");
  
      if (divet==NULL){
        fprintf(stderr, "VH_calls file %s not found.\n", fname); return -1;
      }
  
      if (vcf==NULL){
        fprintf(stderr, "Output file %s cannot be created.\n", outfname); return -1;
      }
  
      if (rep==NULL){
        fprintf(stderr, "Transposons BED file %s not found.\n", rname); return -1;
      }
  
      if (dupfile==NULL){
        fprintf(stderr, "SegDups BED file %s not found.\n", dname); return -1;
      }
  
  
  
      fprintf(vcf, "%s", header);
  
      fprintf(vcf, "##fileDate=%d%s%d%s%d\n", timeinfo->tm_year+1900, (timeinfo->tm_mon+1<10 ? "0" : ""), timeinfo->tm_mon+1, (timeinfo->tm_mday<10 ? "0" : ""), timeinfo->tm_mday);
  
      fprintf(vcf, "%s\n", header2);
      
      for (i=1; i<=22; i++){
  
        sprintf(thischrom, "chr%d", i);
        do_chrom(thischrom, i, minsup, maxsup, mininvsup, maxdup, prune);
      }
  
      i=23;
      sprintf(thischrom, "chrX");
      do_chrom(thischrom, i, minsup, maxsup, mininvsup, maxdup, prune);
  
      if (params->sample_gender == MALE){
      	i=24;
        sprintf(thischrom, "chrY");
        do_chrom(thischrom, i, minsup, maxsup, mininvsup, maxdup, prune);
      }
  
      fprintf(stdout, "Sample = %s\n", samplename);
      fprintf(stdout, "%d dels, %d ins, %d invs.\n", del_cnt, ins_cnt, inv_cnt);
      fprintf(stdout, "Total del %d bp.\nTotal ins %d bp.\n", total_del_bp, total_ins_bp);
      fprintf(stdout, "Min del %d bp.\nMax del %d bp.\n", min_del, max_del);
  
	} 

float calc_dup(int s, int e){
    int len = e-s+1;
    int dupbp = 0;
    int i;
  
      for (i=s-1; i<e; i++) if (dups[i]==1) dupbp++;
  
      if (len > 0)
        return (float)dupbp / (float)len;
      else
      return (dups[s-1]==1) * 100.0;
}

float calc_rep(int s, int e){
    int len = e-s+1;
    int repbp = 0;
    int i;
  
      for (i=s-1; i<e; i++) if (repeats[i]==1) repbp++;
  
      if (len > 0)
        return (float)repbp / (float)len;
      else
      return (repeats[s-1]==1) * 100.0;
  
}


void do_chrom(char *thischrom, int i, int minsup, int maxsup, int mininvsup, float maxdup, float prune){

    int j;
    int s, e;
    char chrom[50];
  
    int os, oe, is, ie;
    char type;
    int support;
    char word[100];
    float dup_content;
    float rep_content_1;
    float rep_content_2;

    char _svtype[20];
    char line[10000];
    float thisprune;

    fprintf(stderr, "chromosome %s\n", thischrom);
    rewind(rep); 
    rewind(dupfile);
    rewind (divet);   
    memset (repeats, 0, MAX*sizeof(char));
    memset (dups, 0, MAX*sizeof(char));

    while (fscanf(rep, "%s\t%d\t%d\n", chrom, &s, &e) > 0){
   		if (!strcmp(chrom, thischrom)){
   			for (j=s;j<e;j++) repeats[j]=1; 
   		}
    }
  
     while (fscanf(dupfile, "%s\t%d\t%d\n", chrom, &s, &e) > 0){
        	if (!strcmp(chrom, thischrom)){
        		for (j=s;j<e;j++) dups[j]=1;
        	} 
    }

    while (fscanf(divet, "%s", word) > 0){
    	if (!strstr(word, "Chr:")){
    		fgets(line, 10000, divet); 
    		continue; 
    	}
	    strcpy (chrom, word+strlen("Chr:"));
	         
	    if (strcmp(chrom, thischrom)) {
	    	fgets(line, 10000, divet); 
	    	continue; 
    	}
     
            
	    fscanf(divet, "%s", word);
	    os = atoi(word + strlen("Start_Outer:"));
	    fscanf(divet, "%s", word);
	    is = atoi(word + strlen("Start_Inner:"));
	    fscanf(divet, "%s", word);
	    ie = atoi(word + strlen("End_Inner:"));
	    fscanf(divet, "%s", word);
	    oe = atoi(word + strlen("End_Outer:"));
	  
		fscanf(divet, "%s", word);
		type = word[7];

		fscanf(divet, "%s", word);
		support = atoi(word + strlen("sup:"));
	  
		if (type=='D' && prune!=10000)
			fscanf(divet, "%s", word); // avg_span

		fscanf(divet, "%s", word); // sum_weight
		fscanf(divet, "%s", word); // lib

		if (!strcmp(samplename, "VHcall_")){
			strcpy(samplename, word+strlen("Lib:"));
		}

		thisprune = 0;
		if (prune == 10000 || type != 'D')
			fgets (word, 100, divet); // libsup; skip
		else{
			fscanf(divet, "%s", word); // libsup
			fscanf(divet, "%s", word); // numcon
			thisprune = atof(word + strlen("NumCon:"));
			//fprintf(stderr, "prune: %f\tthisprune: %f\n", prune, thisprune);
		}

		if (support < minsup) continue;
		if (support > maxsup) continue;
		if (type == 'V' && support < mininvsup) continue;

		if (type == 'D' && thisprune > prune) continue;

		if (type=='D') strcpy(_svtype, "<DEL>");
		else if (type=='I') strcpy(_svtype, "<INS>");
		else if (type=='V') strcpy(_svtype, "<INV>");

		if (type != 'V'){
			dup_content = calc_dup(is, ie);

			if (dup_content > maxdup) continue;	   
		}

		if (type == 'V' || (ram_for_del && type == 'D')){

			rep_content_1 = calc_rep(os, is);
			rep_content_2 = calc_rep(ie, oe);

			if ((rep_content_1 >= 0.75 && rep_content_2 < 0.75) || (rep_content_1 < 0.75 && rep_content_2 >= 0.75)) continue;
			if (type == 'V' && rep_content_1 >= 0.5 && rep_content_2 >= 0.5) continue;

			rep_content_1 = calc_dup(os, is);
			rep_content_2 = calc_dup(ie, oe);

			if ((rep_content_1 >= 0.75 && rep_content_2 < 0.75) || (rep_content_1 < 0.75 && rep_content_2 >= 0.75)) continue;
			if (type == 'V' && rep_content_1 >= 0.5 && rep_content_2 >= 0.5) continue;

		}

		if( (ie-is+1) < 0) continue;

		/* note: the 'N' below should change to the actual basepair divet the reference. See the TODO list at the top of this code */
		fprintf(vcf, "%s\t%d\t%s_%d\t%c\t%s\t.\t.\t", thischrom+strlen("chr"), is, samplename, id++, 'N', _svtype);

		fprintf(vcf, "CIEND=%d,%d;", 0, (oe-ie));
		fprintf(vcf, "CIPOS=%d,%d;", (os-is), 0);
		fprintf(vcf, "END=%d;", ie);
		if (is != ie)
		fprintf(vcf, "IMPRECISE;");

		fprintf(vcf, "SVLEN=%d;", (ie-is+1));
		if (type=='D'){
		fprintf(vcf, "SVTYPE=DEL;");
		total_del_bp += (ie-is+1);
		if (min_del > ie-is+1) min_del = ie-is+1;
		if (max_del < ie-is+1) max_del = ie-is+1;
		del_cnt++;
		}
		else if (type=='I'){
		fprintf(vcf, "SVTYPE=INS;");
		total_ins_bp += (ie-is+1);
		ins_cnt++;
		}
		else if (type=='V'){
		fprintf(vcf, "SVTYPE=INV;");
		inv_cnt++;
		}

		if (strcmp(samplename, "VHcall_"))
		fprintf(vcf, "SAMPLE=%s;", samplename);
		fprintf(vcf, "SVMETHOD=RP;SVALG=VariationHunter;SUP=%d\n", support);
	}
}

int divettovcf( parameters *params, bam_info ** in_bams)
{
	int i, j;
	char cmdline[4096];
	int return_value;

	for( i = 0; i < params->num_bams; i++)
	{
		for ( j = 0; j < in_bams[i]->num_libraries; j++)
	        {
				 conversion(params,in_bams,i,j);
		}
	}  
	
	return RETURN_SUCCESS;
}
