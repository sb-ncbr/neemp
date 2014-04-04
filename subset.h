/*
 * NEEMP - subset.h
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013, 2014
 *
 * */

#ifndef __SUBSET_H__
#define __SUBSET_H__

#include "bitarray.h"

struct stats {

	float R;
	float RMSD;
	float MSE;
	float D_avg;
	float D_max;
};

struct kappa_data {

	const struct subset *parent_subset;

	float kappa;

	float *charges;

	float *parameters_alpha;
	float *parameters_beta;

	struct stats full_stats;
	struct stats *per_at_stats;
	struct stats *per_molecule_stats;
};

void kd_init(struct kappa_data * const kd);
void kd_destroy(struct kappa_data * const kd);
void kd_print_stats(const struct kappa_data * const kd);

float kd_sort_by_return_value(const struct kappa_data * const kd);
int kd_sort_by_is_better(const struct kappa_data * const kd1, const struct kappa_data * const kd2);

struct subset {

	struct bit_array molecules;

	int kappa_data_count;
	struct kappa_data *data;

	/* Pointer to the best kappa */
	struct kappa_data *best;

	struct subset *parent;
};

void ss_destroy(struct subset * const ss);
void print_results(const struct subset * const ss);

#endif /* __SUBSET_H__ */
