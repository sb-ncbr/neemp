/*
 * NEEMP - diffevolution.c
 *
 * by Jana Pazurikova (pazurikova@ics.muni.cz)
 * 2016
 *
 * */

#include <assert.h>
#include <math.h>
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

void generate_random_population(struct subset* ss, float *bounds);
void evolve_kappa(struct kappa_data* trial, struct kappa_data* x, struct kappa_data* a, struct kappa_data *b, float *bounds, double mutation_constant, double recombination_constant);
float get_random_float(float low, float high);
void kd_copy_parameters(struct kappa_data* from, struct kappa_data* to);

void run_diff_evolution(struct subset * const ss)
{
	//create a set of random points in vector space of kappa_data
	//TODO allocate memory for kappa_data array in separate method according to optimization method, should be called also from kappa.c:find_the_best_parameters
	int population_size = 500;//10*(ts.atom_types_count*2+1);
	printf("DE Generating population of size %d\n", population_size);
	ss->kappa_data_count = population_size;
	ss->data = (struct kappa_data*) calloc(ss->kappa_data_count, sizeof(struct kappa_data));
	if(!ss->data)
		EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for kappa data array.\n");
	for (int i = 0; i< ss->kappa_data_count; i++)
	{
		kd_init(&ss->data[i]);
	}
	//set bounds for parameters
	//TODO this should be parameter, ideally possible to set from settings or as described in python prototype by atom type
	float bounds[6]; //[kappa_low, kappa_high, alpha_low, alpha_high, beta_low, beta_high]
	bounds[0] = 0.0005;
	bounds[1] = 3.5;
	bounds[2] = -0.2;
	bounds[3] = 0.5;
	bounds[4] = -0.2;
	bounds[5] = 0.5;
	generate_random_population(ss, bounds);

	//evaluate the fitness function for all points and assign the best
	printf("DE Evaluating fitness function for whole population\n");
	for (int i = 0; i < ss->kappa_data_count; i++)
	{
		calculate_charges(ss, &ss->data[i]);
		calculate_statistics(ss, &ss->data[i]);
	}
	//TODO extract to separate method, also used in kappa.c:find_the_best_parameters
	ss->best = &ss->data[0];
	for (int i = 0; i < ss->kappa_data_count -1; i++)
		if (kd_sort_by_is_better(&ss->data[i], ss->best))
			ss->best = &ss->data[i];

	//run the optimization until converged or for max_iter
	//TODO include iters_max for DE in limits or use one already there (that is used for discard)
	int iter = 0;
	int iters_max = 1000;
	struct kappa_data *trial = (struct kappa_data*)malloc(sizeof(struct kappa_data));
	struct kappa_data *so_far_best = (struct kappa_data*)malloc(sizeof(struct kappa_data));
	kd_init(trial);
	trial->parent_subset = ss;
	kd_init(so_far_best);

	//copy ss->best into so_far_best 
	kd_copy_parameters(ss->best, so_far_best);
	calculate_charges(ss, so_far_best);
	calculate_statistics(ss, so_far_best);

	double mutation_constant = 0.75;
	double recombination_constant = 0.7;
	while (iter < iters_max)
	{
		iter++;
		printf("DE iteration %d\n", iter);
		//select randomly two points from population
		//TODO replace with some real random number generator
		int rand1 = rand() % ss->kappa_data_count;
		int rand2 = rand() % ss->kappa_data_count;
		struct kappa_data a = ss->data[rand1];
		struct kappa_data b = ss->data[rand2];
		mutation_constant = get_random_float(0.5, 1);
		//recombine parts of best, a and b to obtain new trial structure
		evolve_kappa(trial, so_far_best, &a, &b, bounds, mutation_constant, recombination_constant);
		calculate_charges(ss, trial);
		calculate_statistics(ss, trial);
		//if the new structure is better than what we have before, reassign
		if (kd_sort_by_is_better(trial, so_far_best))
		{
			kd_copy_parameters(trial, so_far_best);
			calculate_charges(ss, so_far_best);
			calculate_statistics(ss, so_far_best);
		}
		kd_print_stats(so_far_best);
		print_parameters(so_far_best);
	}
	kd_destroy(trial);
	free(trial);
	kd_copy_parameters(so_far_best, ss->best);
	calculate_charges(ss, ss->best);
	calculate_statistics(ss, ss->best);


}


void generate_random_population(struct subset* ss, float *bounds)
{

	//assign random value to each parameter of each kappa_data
	for (int i = 0; i < ss->kappa_data_count; i++)
	{
		ss->data[i].kappa = get_random_float(bounds[0], bounds[1]);
		for (int j = 0; j < ts.atom_types_count; j++)
		{
			ss->data[i].parameters_alpha[j] = get_random_float(bounds[2], bounds[3]);
			ss->data[i].parameters_beta[j] = get_random_float(bounds[4], bounds[5]);
		}
	}
}

void evolve_kappa(struct kappa_data* trial, struct kappa_data* x, struct kappa_data* a, struct kappa_data *b, float *bounds, double mutation_constant, double recombination_constant)
{
	kd_copy_parameters(x, trial);
	if (get_random_float(0,1) < recombination_constant)
	{
		trial->kappa += mutation_constant*(a->kappa - b->kappa);
		if (bounds[0] > trial->kappa || bounds[1] < trial->kappa)
			trial->kappa = x->kappa;
	}
	for (int i = 0; i < ts.atom_types_count; i++)
		if (get_random_float(0,1) < recombination_constant)
		{
			trial->parameters_alpha[i] += mutation_constant*(a->parameters_alpha[i] - b->parameters_alpha[i]);
			if (bounds[2] > trial->parameters_alpha[i] || bounds[3] < trial->parameters_alpha[i])
				trial->parameters_alpha[i] = x->parameters_alpha[i];
		}
	for (int i = 0; i < ts.atom_types_count; i++)
		if (get_random_float(0,1) < recombination_constant)
		{
			trial->parameters_beta[i] += mutation_constant*(a->parameters_beta[i] - b->parameters_beta[i]);
			if (bounds[4] > trial->parameters_beta[i] || bounds[5] < trial->parameters_beta[i])
				trial->parameters_beta[i] = x->parameters_beta[i];
		}


}

float get_random_float(float low, float high)
{
	//TODO return something closer to random value, this is really pathetic
	float n = low + (float)(rand())/((float)(RAND_MAX/(high-low)));
	return n;
}

void kd_copy_parameters(struct kappa_data* from, struct kappa_data* to)
{
	to->kappa = from->kappa;
	for (int i = 0; i < ts.atom_types_count; i++)
	{
		to->parameters_alpha[i] = from->parameters_alpha[i];
		to->parameters_beta[i] = from->parameters_beta[i];
	}
}
