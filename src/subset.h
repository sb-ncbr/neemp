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
	float R2;
	float R_w;
	float spearman;
	float RMSD;
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
float kd_sort_by_return_value_per_atom(const struct kappa_data * const kd, int i);
int kd_sort_by_is_better(const struct kappa_data * const kd1, const struct kappa_data * const kd2);
void kd_sort_by_is_better_per_atom(int* results_per_atom, const struct kappa_data * const kd1, const struct kappa_data * const kd2, float threshold);

struct subset {

	struct bit_array molecules;

	int kappa_data_count;
	struct kappa_data *data;

	/* Pointer to the best kappa */
	struct kappa_data *best;

	const struct subset *parent;
};

void ss_init(struct subset * const ss, const struct subset * const parent);
void ss_destroy(struct subset * const ss);
void print_results(const struct subset * const ss);
void print_parameters(const struct kappa_data * const kd);
#endif /* __SUBSET_H__ */
