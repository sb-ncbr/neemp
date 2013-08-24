/*
 * NEEMP - subset.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#include <assert.h>
#include <stdlib.h>

#include "bitarray.h"
#include "neemp.h"
#include "settings.h"
#include "structures.h"
#include "subset.h"

extern const struct settings s;
extern const struct training_set ts;

/* Allocate memory for the contents of kappa_data structure */
void kd_init(struct kappa_data * const kd) {

	assert(kd != NULL);

	kd->parameters_alpha = (float *) malloc(ts.atom_types_count * sizeof(float));
	kd->parameters_beta = (float *) malloc(ts.atom_types_count * sizeof(float));
	if(!kd->parameters_alpha || !kd->parameters_beta)
		EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for parameters array.\n");

	kd->charges = (float *) malloc(ts.atoms_count * sizeof(float));
	if(!kd->charges)
		EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for charges array.\n");

	kd->max_D_per_atom_type = (float *) calloc(ts.atom_types_count, sizeof(float));
	kd->avg_D_per_atom_type = (float *) calloc(ts.atom_types_count, sizeof(float));

	if(!kd->max_D_per_atom_type || !kd->avg_D_per_atom_type)
		EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for stats array.\n");
}

/* Destroy contents of the kappa_data structure */
void kd_destroy(struct kappa_data * const kd) {

	assert(kd != NULL);

	free(kd->parameters_alpha);
	free(kd->parameters_beta);
	free(kd->charges);
	free(kd->max_D_per_atom_type);
	free(kd->avg_D_per_atom_type);
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

	printf("\nResults\n\n");

	printf("Used molecules: %5d\n", b_count_bits(&ss->molecules));
	printf("K: %6.4f | R: %6.4f   RMSD: %6.4f   MSE: %6.4f   D_avg: %6.4f   D_max: %6.4f\n",
		ss->best->kappa, ss->best->R, ss->best->RMSD, ss->best->MSE, ss->best->D_avg, ss->best->D_max);

	printf("Atom type             A       B\t\t  D_max\t  D_avg\n");
	for(int i = 0; i < ts.atom_types_count; i++) {
		printf(" %2s %1d   \t%7.4f\t%7.4f\t\t%7.3f\t%7.3f\n",
			convert_Z_to_symbol(ts.atom_types[i].Z), ts.atom_types[i].bond_order,
			ss->best->parameters_alpha[i], ss->best->parameters_beta[i],
			ss->best->max_D_per_atom_type[i], ss->best->avg_D_per_atom_type[i]);
	}

	printf("\n");
}

/* Return the value which we selected for sorting */
float kd_sort_by_return_value(const struct kappa_data * const kd) {

	switch(s.sort_by) {
		case SORT_R:
			return kd->R;
		case SORT_RMSD:
			return kd->RMSD;
		case SORT_MSE:
			return kd->MSE;
		case SORT_D_AVG:
			return kd->D_avg;
		case SORT_D_MAX:
			return kd->D_max;
		default:
			/* Something bad happened */
			assert(0);
	}
}

/* Determine if kd1 is better than kd2 in terms of the sort-by value */
int kd_sort_by_is_better(const struct kappa_data * const kd1, const struct kappa_data * const kd2) {

	if(s.sort_by == SORT_R)
		return kd1->R > kd2->R;
	else
		return kd_sort_by_return_value(kd1) < kd_sort_by_return_value(kd2);
}
