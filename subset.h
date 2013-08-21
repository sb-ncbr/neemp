/*
 * NEEMP - subset.h
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#ifndef __SUBSET_H__
#define __SUBSET_H__

#include "bitarray.h"

struct kappa_data;
struct subset;

struct subset {

	struct bit_array molecules;

	int kappa_data_count;
	struct kappa_data *data;

	/* Pointer to the best kappa */
	struct kappa_data *best;
};

void ss_destroy(struct subset * const ss);

struct kappa_data {

	const struct subset *parent_subset;

	float kappa;

	float *charges;

	float *parameters_alpha;
	float *parameters_beta;

	/* Some statistics */
	float R;			/* Pearson's correlation squared for each molecule, then averaged */
	float RMSD;			/* Root-mean-square deviation for each molecule, then averaged */
	float D;			/* Absolute difference for all atoms */
};

void kd_init(struct kappa_data * const kd);
void kd_destroy(struct kappa_data * const kd);

#endif /* __SUBSET_H__ */
