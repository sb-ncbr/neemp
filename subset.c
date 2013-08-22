/*
 * NEEMP - subset.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#include <assert.h>
#include <stdlib.h>

#include "neemp.h"
#include "structures.h"
#include "subset.h"

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
}
