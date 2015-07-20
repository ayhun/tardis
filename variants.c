#include <stdlib.h>
#include <string.h>

#include "tardis.h"
#include "variants.h"

#define add_to_beginning 1

enum SVTYPE getEnum(char c) {
    switch (c) {
        case 'D':
            return DEL;
        case 'V':
            return INV;
        case 'S':
            return INS;
    }

    return DEL;
}

char * svtypeToChar(enum SVTYPE svt) {
    switch (svt) {
        case DEL:
            return "DEL";
        case INS:
            return "INS";
        case INV:
            return "INV";
    }

}

// create a new structure and return address

struct strvar* new_strvar(int outer_start, int inner_start, int outer_end, int inner_end, enum SVTYPE svtype, int rp_sup, int sr_sup, float avg_edit, int min_svlen, int max_svlen, char *samples, double conf_score) {
    struct strvar* a_strvar = malloc(sizeof (struct strvar));

    a_strvar->outer_start = outer_start;
    a_strvar->inner_start = inner_start;
    a_strvar->outer_end = outer_end;
    a_strvar->inner_end = inner_end;
    a_strvar->svtype = svtype;
    a_strvar->rp_sup = rp_sup;
    a_strvar->sr_sup = sr_sup; //******
    a_strvar->avg_edit = avg_edit;
    a_strvar->min_svlen = min_svlen;
    a_strvar->max_svlen = max_svlen;
    a_strvar->samples = samples; //******
    a_strvar->conf_score = conf_score; //******
    a_strvar->next = NULL;
    a_strvar->head = NULL;

    return a_strvar;
}

struct strvar ** init_vars(int num_chroms) {
    struct strvar** variations;
    // use calloc to ensure null values at memory allocation
    variations = (struct strvar**) calloc(num_chroms+1, sizeof (struct strvar*));

    return variations;
}

int get_index(struct strvar ** variations, char* chr, int create_index_if_not_exist){
    if(variations[0]==NULL){
        struct chr_index* a_chr_index = malloc(sizeof (struct chr_index));
        strcpy(a_chr_index->chrom_name, chr);
        a_chr_index->index = 1;
        a_chr_index->next = NULL;
        variations[0] = (struct strvar*) a_chr_index;
        return a_chr_index->index;
    } else {
        struct chr_index* curr = *((struct chr_index**) variations);
        while( (strcmp(curr->chrom_name, chr) != 0) && curr->next != NULL)
            curr = curr->next;

        if(strcmp(curr->chrom_name, chr) == 0){ //while ended because we found chrom_name
            return curr->index;
        } else if (curr->next == NULL && create_index_if_not_exist == 1){ //while ended because list ended, and we can create new index
            curr->next = malloc(sizeof (struct chr_index));
            strcpy(curr->next->chrom_name, chr);
            curr->next->index = curr->index + 1;
            curr->next->next = NULL;
            return curr->next->index;
        } else {
            return -100;//we can't create new index
        }

    }

}

void add_strvar(struct strvar ** variations, char* chr_, struct strvar* sv) {
    int chr = get_index(variations, chr_,1);
    //struct strvar* ptr = variations[chr];
    struct strvar* ptr = (struct strvar*) (variations + chr);

    if (variations[chr] == NULL) {
        sv->next = NULL; //ensure last node doesn't have the next pointer
        sv->head = variations + chr;
        variations[chr] = sv;
        return;
    }

    if (add_to_beginning) {
        sv->next = variations[chr];
        sv->head = variations + chr;
        variations[chr] = sv;
    } else {
        while (ptr->next != NULL)
            ptr = ptr->next;

        ptr->next = sv;
        sv->next = NULL;
        sv->head = variations + chr;
    }
}

void print_all(struct strvar ** variations) {
    struct chr_index* curr = *((struct chr_index**) variations);
    while(curr != NULL){
        print_chr(variations, (char*) &(curr->chrom_name));
        curr = curr->next;
    }
}

void print_chr(struct strvar ** variations, char* chr_) {
    int chr = get_index(variations, chr_,0);

    if (chr < 1 || variations[chr] == NULL) {
        return;
    }

    struct strvar* curr_strvar = variations[chr];

    while (curr_strvar != NULL) {
        printf("Chr:%-11s Start_Outer:%-11i Start_Inner:%-11i End_Outer:%-11i End_Inner:%-11i SVtype:%-3s rp_sup:%-3i sr_sup:%-1i min svlen:%-6i max svlen:%-6i confidence:%f\n",
                chr_, curr_strvar->outer_start, curr_strvar->inner_start,
                curr_strvar->outer_end, curr_strvar->inner_end, svtypeToChar(curr_strvar->svtype), curr_strvar->rp_sup, curr_strvar->sr_sup,
                curr_strvar->min_svlen, curr_strvar->max_svlen, curr_strvar->conf_score);
        curr_strvar = curr_strvar->next;
    }
}

void free_all(struct strvar ** variations) {
    struct chr_index* curr = *((struct chr_index**) variations);
    struct chr_index* tmp;

    while(curr != NULL){
        free_chr(variations, (char*) &(curr->chrom_name));
        tmp = curr;
        curr = curr->next;
        free(tmp);
    }

    free(variations);
}

void free_chr(struct strvar ** variations, char* chr_) {
    int chr = get_index(variations, chr_, 0);

    if (chr < 1 || variations[chr] == NULL) {
        return;
    }

    struct strvar* head = variations[chr];
    struct strvar* curr;
    while ((curr = head) != NULL) {
        head = head->next;
        free(curr);
    }
}

void remove_strvar(struct strvar * variation) {
    struct strvar* prev = *(variation->head);

    // first element is the element to delete
    if(prev == variation){
        struct strvar** tmp = variation->head;
        *tmp = variation->next;
        free(variation);
    } else {
        while (prev != NULL && prev->next != variation)
            prev = prev->next;

        if(prev != NULL){
            prev->next = variation->next;
            free(variation);
        }
    }
}
