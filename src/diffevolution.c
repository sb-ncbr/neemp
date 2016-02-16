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
#include "latin_random.h"
#include "diffevolution.h"

extern const struct training_set ts;
extern const struct settings s;

void generate_random_population(struct subset* ss, float *bounds);
int evolve_kappa(struct kappa_data* trial, struct kappa_data* x, struct kappa_data* a, struct kappa_data *b, float *bounds, double mutation_constant, double recombination_constant);
void compute_parameters_bounds(float* bounds, int by_atom_type);
float get_random_float(float low, float high);
float interpolate_to_different_bounds(float x, float low, float high);
void kd_copy_parameters(struct kappa_data* from, struct kappa_data* to);

void run_diff_evolution(struct subset * const ss)
{
	//create a set of random points in vector space of kappa_data
	//TODO allocate memory for kappa_data array in separate method according to optimization method, should be called also from kappa.c:find_the_best_parameters
	printf("DE Generating population of size %d\n", s.population_size);
	ss->kappa_data_count = s.population_size;
	ss->data = (struct kappa_data*) calloc(ss->kappa_data_count, sizeof(struct kappa_data));
	if(!ss->data)
		EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for kappa data array.\n");
	for (int i = 0; i< ss->kappa_data_count; i++)
	{
		kd_init(&ss->data[i]);
	}
	//set bounds for parameters
	float *bounds = (float*) malloc((ts.atom_types_count*2+1)*2*sizeof(float));
	compute_parameters_bounds(bounds, 0);
	generate_random_population(ss, bounds);

	//evaluate the fitness function for all points and assign the best
	printf("DE Evaluating fitness function for whole population\n");
	for (int i = 0; i < ss->kappa_data_count; i++)
	{
		calculate_charges(ss, &ss->data[i]);
		calculate_statistics_by_sort_mode(ss, &ss->data[i]);
	}
	//TODO extract to separate method, also used in kappa.c:find_the_best_parameters
	ss->best = &ss->data[0];
	for (int i = 0; i < ss->kappa_data_count -1; i++)
		if (kd_sort_by_is_better(&ss->data[i], ss->best))
			ss->best = &ss->data[i];

	//run the optimization until converged or for max_iter
	//TODO include iters_max for DE in limits or use one already there (that is used for discard)
	int iter = 0;
	struct kappa_data *trial = (struct kappa_data*)malloc(sizeof(struct kappa_data));
	struct kappa_data *so_far_best = (struct kappa_data*)malloc(sizeof(struct kappa_data));
	kd_init(trial);
	trial->parent_subset = ss;
	kd_init(so_far_best);

	//copy ss->best into so_far_best 
	kd_copy_parameters(ss->best, so_far_best);
	calculate_charges(ss, so_far_best);
	calculate_statistics(ss, so_far_best);
	float mutation_constant = s.mutation_constant;
	int iters_with_evolution=0;
	while (iter < s.limit_de_iters)
	{
		iter++;
		printf("DE iteration %d\n", iter);
		//select randomly two points from population
		//TODO replace with some real random number generator
		int rand1 = rand() % ss->kappa_data_count;
		int rand2 = rand() % ss->kappa_data_count;
		struct kappa_data a = ss->data[rand1];
		struct kappa_data b = ss->data[rand2];
		if (s.dither)
			mutation_constant = get_random_float(0.5, 1);
		//recombine parts of best, a and b to obtain new trial structure
		if (evolve_kappa(trial, so_far_best, &a, &b, bounds, mutation_constant, s.recombination_constant))
		{
			iters_with_evolution++;
			calculate_charges(ss, trial);
			calculate_statistics_by_sort_mode(ss, trial);
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
	}
	kd_destroy(trial);
	free(trial);
	free(bounds);
	kd_copy_parameters(so_far_best, ss->best);
	calculate_charges(ss, ss->best);
	calculate_statistics(ss, ss->best);
	printf("From %d iterations, %d introduced a change\n", s.limit_de_iters, iters_with_evolution);


}


void generate_random_population(struct subset* ss, float *bounds)
{
	//get random number by Latin Hypercube Sampling
	int dimensions_count = ts.atom_types_count*2+1;
	int points_count = s.population_size;
	int seed = 100;
	double* random_lhs = latin_random_new(dimensions_count, points_count, &seed);

	//redistribute random_lhs[dim_num, point_num] to ss->data
	for (int i = 0; i < points_count; i++) {
		ss->data[i].kappa = interpolate_to_different_bounds(random_lhs[i*dimensions_count], bounds[0], bounds[1]);
		for (int j = 0; j < ts.atom_types_count; j++) {
			//interpolate random numbers from [0,1] to new bounds
			ss->data[i].parameters_alpha[j] = interpolate_to_different_bounds(random_lhs[i*dimensions_count + 1 + j*2], bounds[2 + j*4], bounds[2 + j*4 + 1]);
			ss->data[i].parameters_beta[j] = interpolate_to_different_bounds(random_lhs[i*dimensions_count + 2 + j*2 ], bounds[2 + j*4 + 2], bounds[2 + j*4 + 3]);
		}
	}

}

int evolve_kappa(struct kappa_data* trial, struct kappa_data* x, struct kappa_data* a, struct kappa_data *b, float *bounds, double mutation_constant, double recombination_constant)
{
	int changed = 0;
	kd_copy_parameters(x, trial);
	for (int i = 0; i < ts.atom_types_count; i++)
		if (get_random_float(0,1) < recombination_constant)
		{
			changed++;
			trial->parameters_alpha[i] += mutation_constant*(a->parameters_alpha[i] - b->parameters_alpha[i]);
			if (bounds[2 + i*4] > trial->parameters_alpha[i] || bounds[2 + i*4 + 1] < trial->parameters_alpha[i]) {
				trial->parameters_alpha[i] = x->parameters_alpha[i];
				changed--;
			}
		}
	for (int i = 0; i < ts.atom_types_count; i++)
		if (get_random_float(0,1) < recombination_constant)
		{
			changed++;
			trial->parameters_beta[i] += mutation_constant*(a->parameters_beta[i] - b->parameters_beta[i]);
			if (bounds[2 + i*4 + 2] > trial->parameters_beta[i] || bounds[2 + i*4 + 3] < trial->parameters_beta[i]) {
				trial->parameters_beta[i] = x->parameters_beta[i];
				changed--;
			}
		}

	//always change kappa
	trial->kappa += mutation_constant*(a->kappa - b->kappa)/(bounds[1]-bounds[0]);
	changed++;
	if (bounds[0] > trial->kappa || bounds[1] < trial->kappa) {
		trial->kappa -= mutation_constant*(a->kappa - b->kappa)/(bounds[1]-bounds[0]);
	}
	
	return changed;
}

void compute_parameters_bounds(float* bounds, int by_atom_type) {          
	//returns bounds[kappa_low, kappa_high, alpha_1_low, alpha_1_high, beta_1_low, beta_1_high, alpha_2_low, ...]
	float toH = 0.036749309;
	bounds[0] = 0.0005; //kappa_low
	bounds[1] = 3.5; //kappa_high
	for (int j = 0; j < ts.atom_types_count; j++) {	
		if (by_atom_type) {
			//set bounds for particular atom type
			int atom_number = ts.atom_types[j].Z;
			bounds[2 + j*4] = (ionenergies[atom_number] + affinities[atom_number])/2*toH - 0.1; //alpha_low
			bounds[2 + j*4 + 1] = bounds[2 + j*4] + 0.2; //alpha_high
			bounds[2 + j*4 + 2] = (ionenergies[atom_number] - affinities[atom_number])/2*toH - 0.1; //beta_low
			bounds[2 + j*4 + 3]	= bounds[2 + j*4 + 2] + 0.2; //beta_high
		}
		else {
			bounds[2 + j*4] = 2;
			bounds[2 + j*4 + 1] = 3;
			bounds[2 + j*4 + 2] = 0;
			bounds[2 + j*4 + 3] = 0.8;
		}
	}
}

float get_random_float(float low, float high) {
	float n = low + (float)(rand())/((float)(RAND_MAX/(high-low)));
	return n;
}

float interpolate_to_different_bounds(float x, float low, float high) {
	return low + x*(high-low);
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
