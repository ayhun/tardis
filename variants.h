#ifndef __VARIANTS
#define __VARIANTS

#include <stdio.h>

enum SVTYPE {
    DEL, INS, INV, TANDUP, MEI, TRANSCHROM
};
enum SVTYPE getEnum(char c);
char * svtypeToChar(enum SVTYPE svt);

typedef struct strvar {
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
    struct strvar **head; /* head pointer points to pointer that points to head element */

}strvar;

typedef struct chr_index {
    char chrom_name[20];
    int index;
    struct chr_index *next;
}chr_index;
/* functions */

struct strvar ** init_vars(int num_chroms);
struct strvar* new_strvar(int outer_start, int inner_start, int outer_end, int inner_end, enum SVTYPE svtype, int rp_sup, int sr_sup, float avg_edit, int min_svlen, int max_svlen, char *samples, double conf_score);
void add_strvar(struct strvar ** variations, char* chr, struct strvar* sv);
void remove_strvar(struct strvar * variation);
void print_all(struct strvar ** variations);
void print_chr(struct strvar ** variations, char* chr);
void free_all(struct strvar ** variations);
void free_chr(struct strvar ** variations, char* chr);
#endif
