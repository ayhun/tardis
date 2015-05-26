#ifndef __VARIANTS
#define __VARIANTS

enum SVTYPE {DEL, INS, INV, TANDUP, MEI, TRANSCHROM};

struct strvar
{
  int outer_start; /* outer start coordinate */
  int inner_start; /* inner start coordinate */
  int outer_end; /* outer end coordinate */
  int inner_end; /* inner end coordinate */
  enum SVTYPE svtype; /* type of SV */
  int rp_sup; /* number of RP support */
  int sr_sup; /* number of SR support */
  float avg_edit; /* average edit distance of reads mapped to reference */
  int min_svlen; /* lower bound of SV size */
  int max_svlen; /* upper bound of SV size */
  char *samples; /* list of samples that carry the SV. Might need a separate linked list */
  double conf_score; /* confidence score for the called variant */
   
  struct strvar *next; /* next pointer for linked list */

};

/* functions */

/*
struct strvar ** init_vars(int num_chroms);
int  add_strvar
int  remove_strvar
int  edit_strvar
 */

#endif
