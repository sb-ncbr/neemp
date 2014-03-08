/*
 * NEEMP - subset.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013, 2014
 *
 * */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "bitarray.h"
#include "neemp.h"
#include "settings.h"
#include "structures.h"
#include "subset.h"

extern const struct settings s;
extern const struct training_set ts;


/* Allocate memory for the contents of statistics structure */
void st_init(struct statistics * const st) {

	assert(st != NULL);

	st->max_D_per_atom_type = (float *) calloc(ts.atom_types_count, sizeof(float));
	st->avg_D_per_atom_type = (float *) calloc(ts.atom_types_count, sizeof(float));
	st->R_per_atom_type = (float *) calloc(ts.atom_types_count, sizeof(float));
	st->RMSD_per_atom_type = (float *) calloc(ts.atom_types_count, sizeof(float));
	st->MSE_per_atom_type = (float *) calloc(ts.atom_types_count, sizeof(float));

	if(!st->max_D_per_atom_type || !st->avg_D_per_atom_type || !st->R_per_atom_type || \
		!st->RMSD_per_atom_type || !st->MSE_per_atom_type)
		EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for statistics.\n");
}

/* Destroy contents of the statistics structure */
void st_destroy(struct statistics * const st) {

	assert(st != NULL);

	free(st->R_per_atom_type);
	free(st->RMSD_per_atom_type);
	free(st->MSE_per_atom_type);
	free(st->max_D_per_atom_type);
	free(st->avg_D_per_atom_type);
}

/* Allocate memory for the contents of kappa_data structure */
void kd_init(struct kappa_data * const kd) {

	assert(kd != NULL);

	kd->parameters_alpha = (float *) calloc(ts.atom_types_count, sizeof(float));
	kd->parameters_beta = (float *) calloc(ts.atom_types_count, sizeof(float));
	if(!kd->parameters_alpha || !kd->parameters_beta)
		EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for parameters array.\n");

	kd->charges = (float *) malloc(ts.atoms_count * sizeof(float));
	if(!kd->charges)
		EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for charges array.\n");

	st_init(&kd->stats);
}

/* Destroy contents of the kappa_data structure */
void kd_destroy(struct kappa_data * const kd) {

	assert(kd != NULL);

	st_destroy(&kd->stats);

	free(kd->parameters_alpha);
	free(kd->parameters_beta);
	free(kd->charges);
}

/* Destroy contents of the subset */
void ss_destroy(struct subset * const ss) {

	assert(ss != NULL);

	for(int i = 0; i < ss->kappa_data_count; i++)
		kd_destroy(&ss->data[i]);

	b_destroy(&ss->molecules);
	free(ss->data);
}

/* Prints the parameters and associated stats */
void print_results(const struct subset * const ss) {

	assert(ss != NULL);
	assert(ss->best != NULL);

	printf("Used molecules: %5d\n", b_count_bits(&ss->molecules));
	kd_print_stats(ss->best);

	printf("Atom type             A       B\t\t     R\t    RMSD       MSE\t  D_avg\t  D_max\n");
	for(int i = 0; i < ts.atom_types_count; i++) {
		char buff[10];
		at_format_text(&ts.atom_types[i], buff);
		printf(" %s   \t%7.4f\t%7.4f\t\t%6.4f  %4.2e  %4.2e\t%7.3f\t%7.3f\n",
			buff, ss->best->parameters_alpha[i], ss->best->parameters_beta[i],
			ss->best->stats.R_per_atom_type[i], ss->best->stats.RMSD_per_atom_type[i], ss->best->stats.MSE_per_atom_type[i],
			ss->best->stats.avg_D_per_atom_type[i], ss->best->stats.max_D_per_atom_type[i]);
	}

	printf("\n");
}

/* Return the value which we selected for sorting */
float kd_sort_by_return_value(const struct kappa_data * const kd) {

	switch(s.sort_by) {
		case SORT_R:
			return kd->stats.R;
		case SORT_RMSD:
			return kd->stats.RMSD;
		case SORT_MSE:
			return kd->stats.MSE;
		case SORT_D_AVG:
			return kd->stats.D_avg;
		case SORT_D_MAX:
			return kd->stats.D_max;
		default:
			/* Something bad happened */
			assert(0);
	}
}

/* Determine if kd1 is better than kd2 in terms of the sort-by value */
int kd_sort_by_is_better(const struct kappa_data * const kd1, const struct kappa_data * const kd2) {

	if(s.sort_by == SORT_R)
		return kd1->stats.R > kd2->stats.R;
	else
		return kd_sort_by_return_value(kd1) < kd_sort_by_return_value(kd2);
}

/* Print all the statistics for the particular kappa data */
void kd_print_stats(const struct kappa_data * const kd) {

	/* Allocate buffer large enough to hold the entire message */
	char message[100];
	memset(message, 0, 100 * sizeof(char));

	snprintf(message, 100, "K: %6.4f |  R: %6.4f  RMSD: %4.2e  MSE: %4.2e  D_avg: %4.2e  D_max: %4.2e\n",
		kd->kappa, kd->stats.R, kd->stats.RMSD, kd->stats.MSE, kd->stats.D_avg, kd->stats.D_max);

	/* Print '*' at the position of the stat which is used for sorting */
	switch(s.sort_by) {
		case SORT_R:
			message[12] = '*';
			break;
		case SORT_RMSD:
			message[23] = '*';
			break;
		case SORT_MSE:
			message[39] = '*';
			break;
		case SORT_D_AVG:
			message[54] = '*';
			break;
		case SORT_D_MAX:
			message[71] = '*';
			break;
		default:
			/* Something bad happened */
			assert(0);
	}

	printf("%s", message);
}
