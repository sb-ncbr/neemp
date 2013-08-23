/*
 * NEEMP - kappa.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "eem.h"
#include "kappa.h"
#include "neemp.h"
#include "parameters.h"
#include "settings.h"
#include "subset.h"
#include "statistics.h"
#include "structures.h"

extern const struct training_set ts;
extern const struct settings s;

static void full_scan(struct subset * const ss);
static void brent(struct subset * const ss);
static void perform_calculations(struct subset * const ss, struct kappa_data * const kd);

/* Perform all three step for one value of kappa */
static void perform_calculations(struct subset * const ss, struct kappa_data * const kd) {

	assert(ss != NULL);
	assert(kd != NULL);

	calculate_parameters(ss, kd);
	calculate_charges(ss, kd);
	calculate_statistics(ss, kd);
}

static void full_scan(struct subset * const ss) {

	assert(ss != NULL);

	const int steps = s.full_scan_only ? ss->kappa_data_count : ss->kappa_data_count - 1;

	for(int i = 0; i < steps && i * s.full_scan_precision < s.kappa_max; i++) {
		ss->data[i].kappa = i * s.full_scan_precision;
		perform_calculations(ss, &ss->data[i]);
		printf("%6.4f | R: %8.4f   RMSD: %8.4f   MSE: %9.4f   D: %8.4f\n",
			ss->data[i].kappa, ss->data[i].R, ss->data[i].RMSD, ss->data[i].MSE, ss->data[i].D);
	}
}

static void brent(struct subset * const ss) {

	assert(ss != NULL);

	full_scan(ss);
}

/* Finds the best parameters for a particular subset */
void find_the_best_parameters_for_subset(struct subset * const ss) {

	assert(ss != NULL);

	if(s.kappa_set > 1e-10) {

		ss->kappa_data_count = 1;
		ss->data = (struct kappa_data *) malloc(1 * sizeof(struct kappa_data));
		kd_init(&ss->data[0]);
		ss->data[0].kappa = s.kappa_set;

		perform_calculations(ss, &ss->data[0]);
	}
	else {
		ss->kappa_data_count = 2 + (int) (s.kappa_max / s.full_scan_precision);
		ss->data = (struct kappa_data *) calloc(ss->kappa_data_count, sizeof(struct kappa_data));
		if(!ss->data)
			EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for kappa data array.\n");

		for(int i = 0; i < ss->kappa_data_count; i++)
			kd_init(&ss->data[i]);

		if(s.full_scan_only)
			full_scan(ss);
		else
			brent(ss);
	}
	float best_R = 0.0f;

	for(int i = 0; i < ss->kappa_data_count; i++) {
		if(ss->data[i].R > best_R) {
			ss->best = &ss->data[i];
			best_R = ss->data[i].R;
		}
	}
}
