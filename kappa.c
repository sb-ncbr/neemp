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

	int steps = 10;
	float kappa_step = 0.1f;

	for(int i = 0; i < steps; i++) {
		ss->data[i].kappa = i * kappa_step;
		perform_calculations(ss, &ss->data[i]);
		printf("K = %6.4f R = %6.4f\n", ss->data[i].kappa, ss->data[i].R);
	}
}

static void brent(struct subset * const ss) {

	assert(ss != NULL);

	full_scan(ss);
}

/* Finds the best parameters for a particular subset */
void find_the_best_parameters_for_subset(struct subset * const ss) {

	assert(ss != NULL);

	if(s.full_scan) {
		ss->kappa_data_count = 10;
		ss->data = (struct kappa_data *) malloc(ss->kappa_data_count * sizeof(struct kappa_data));
		for(int i = 0; i < ss->kappa_data_count; i++)
			kd_init(&ss->data[i]);

		full_scan(ss);
	}
	else {
		ss->kappa_data_count = 10;
		ss->data = (struct kappa_data *) malloc(ss->kappa_data_count * sizeof(struct kappa_data));
		for(int i = 0; i < ss->kappa_data_count; i++)
			kd_init(&ss->data[i]);

		brent(ss);
	}

	float best_R = 0.0f;

	for(int i = 0; i < ss->kappa_data_count; i++) {
		if(ss->data[i].R > best_R) {
			ss->best = &ss->data[i];
			best_R = ss->data[i].R;
		}
	}

//	printf("K = %6.4f R = %6.4f D = %6.4f RMSD = %6.4f\n", ss->best->kappa, ss->best->R, ss->best->D, ss->best->RMSD);
}
