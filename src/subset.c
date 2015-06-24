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

/* Initialize new subset structure from parent */
void ss_init(struct subset * const ss, const struct subset * const parent) {

	ss->weights = (float *) calloc(ts.atom_types_count, sizeof(float));
	if(!ss->weights)
		EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for weights array.\n");

	b_init(&ss->molecules, ts.molecules_count);
	ss->parent = parent;
	if(parent) {
		b_set_as(&ss->molecules, &parent->molecules);
		memcpy(ss->weights, parent->weights, ts.atom_types_count * sizeof(float));
	}
	else
		b_set_all(&ss->molecules);
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

	kd->per_at_stats = (struct stats *) calloc(ts.atom_types_count, sizeof(struct stats));
	kd->per_molecule_stats = (struct stats *) calloc(ts.molecules_count, sizeof(struct stats));
}

/* Destroy contents of the kappa_data structure */
void kd_destroy(struct kappa_data * const kd) {

	assert(kd != NULL);

	free(kd->parameters_alpha);
	free(kd->parameters_beta);
	free(kd->charges);
	free(kd->per_at_stats);
	free(kd->per_molecule_stats);
}

/* Destroy contents of the subset */
void ss_destroy(struct subset * const ss) {

	assert(ss != NULL);

	for(int i = 0; i < ss->kappa_data_count; i++)
		kd_destroy(&ss->data[i]);

	b_destroy(&ss->molecules);
	free(ss->data);
	free(ss->weights);
}

/* Print loaded parameters */
void print_parameters(const struct kappa_data * const kd) {

	assert(kd != NULL);

	printf("\nLoaded parameters:\n");
	printf("K: %6.4f\n", kd->kappa);
	printf("Atom type             A       B\n");
	for(int i = 0; i < ts.atom_types_count; i++) {
		char buff[10];
		at_format_text(&ts.atom_types[i], buff);
		printf(" %s   \t%7.4f\t%7.4f\n", buff, kd->parameters_alpha[i], kd->parameters_beta[i]);
	}

	printf("\n");
}

/* Prints the parameters and associated stats */
void print_results(const struct subset * const ss) {

	assert(ss != NULL);
	assert(ss->best != NULL);

	printf("\nUsed molecules: %5d\n", b_count_bits(&ss->molecules));
	kd_print_stats(ss->best);

	printf("Atom type            A       B           R    RMSD   D_avg   D_max\n");
	for(int i = 0; i < ts.atom_types_count; i++) {
		char buff[10];
		at_format_text(&ts.atom_types[i], buff);
		printf(" %s  \t%6.4f\t%6.4f      %6.4f  %6.4f  %6.4f  %6.4f\n",
			buff, ss->best->parameters_alpha[i], ss->best->parameters_beta[i],
			ss->best->per_at_stats[i].R, ss->best->per_at_stats[i].RMSD,
			ss->best->per_at_stats[i].D_avg, ss->best->per_at_stats[i].D_max);
	}

	printf("\n");
}

/* Return the value which we selected for sorting */
float kd_sort_by_return_value(const struct kappa_data * const kd) {

	switch(s.sort_by) {
		case SORT_R:
			return kd->full_stats.R;
		case SORT_R_WEIGHTED:
			return kd->full_stats.R_weighted;
		case SORT_RMSD:
			return kd->full_stats.RMSD;
		case SORT_D_AVG:
			return kd->full_stats.D_avg;
		case SORT_D_MAX:
			return kd->full_stats.D_max;
		default:
			/* Something bad happened */
			assert(0);
	}
}

/* Determine if kd1 is better than kd2 in terms of the sort-by value */
int kd_sort_by_is_better(const struct kappa_data * const kd1, const struct kappa_data * const kd2) {

	/* Higher R/R_weighted is better, otherwise prefer lower value */
	if(s.sort_by == SORT_R || s.sort_by == SORT_R_WEIGHTED)
		return kd_sort_by_return_value(kd1) > kd_sort_by_return_value(kd2);
	else
		return kd_sort_by_return_value(kd1) < kd_sort_by_return_value(kd2);
}

/* Print all the statistics for the particular kappa data */
void kd_print_stats(const struct kappa_data * const kd) {

	/* Allocate buffer large enough to hold the entire message */
	char message[100];
	memset(message, 0, 100 * sizeof(char));

	snprintf(message, 100, "K: %6.4f |  R: %6.4f  R_w: %6.4f  RMSD: %6.4f  D_avg: %6.4f  D_max: %6.4f\n",
		kd->kappa, kd->full_stats.R, kd->full_stats.R_weighted, kd->full_stats.RMSD, kd->full_stats.D_avg, kd->full_stats.D_max);

	printf("%s", message);
}


/* Print assigned weights for R_w calculation */
void print_weights(const struct subset * const ss) {

	assert(ss != NULL);

	printf("The weights used to calculate R_w:\n");
	for(int i = 0; i < ts.atom_types_count; i++) {
			char buff[10];
			at_format_text(&ts.atom_types[i], buff);
			printf("%s: %8.4f\n", buff, ss->weights[i]);
		}

	printf("\n");
}
