/*
 * NEEMP - guidedmin.c
 *
 * by Jana Pazurikova (pazurikova@ics.muni.cz)
 * 2016
 *
 * */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "omp.h"
#include "eem.h"
#include "kappa.h"
#include "neemp.h"
#include "parameters.h"
#include "settings.h"
#include "subset.h"
#include "statistics.h"
#include "structures.h"
#include "../externals/lhs/latin_random.h"
#include "guidedmin.h"

extern const struct training_set ts;
extern const struct settings s;
extern const float ionenergies[];
extern const float affinities[];


/* Run guided minimization algorithm to find the best set of parameters for calculation of partial charges. */ 
void run_guided_min(struct subset * const ss) {

	/* Create a set of random points in vector space of kappa_data */            
	if (s.verbosity >= VERBOSE_KAPPA)
		printf("GM Generating %d vectors\n", s.gm_size);

	/* Set bounds for each parameters in kappa_data */
	float *bounds = (float*) malloc((ts.atom_types_count*2+1)*2*sizeof(float));
	//compute bounds, 0 means set them to fixed numbers taken from Tomas's full scan, 1 means try to find them with broad search
	compute_parameters_bounds(bounds, 0);
	//generate population
	fill_ss(ss, s.gm_size); 
	generate_random_population(ss, bounds, s.gm_size);
    de_ss = ss;
	//evaluate the fitness function for all points
	if (s.verbosity >= VERBOSE_KAPPA)
		printf("GM Calculating charges and evaluating fitness function for whole set\n");
	int i = 0;
#pragma omp parallel for num_threads(s.gm_threads) default(shared) private(i)
	for (i = 0; i < ss->kappa_data_count; i++) {
		calculate_charges(ss, &ss->data[i]);
		calculate_statistics(ss, &ss->data[i]);
	}

	//minimize part of population that has R>0.3
	int minimized_initial = 0;
	if (s.verbosity >= VERBOSE_KAPPA)
		printf("GM minimizing population\n");
	minimized_initial = minimize_part_of_gm_set(ss, s.gm_iterations_beg);
	
	//if we minimized zero data, inform user and suggest larger set
	if (minimized_initial == 0) {
		EXIT_ERROR(RUN_ERROR, "No vector in set was worth minimizing. Please choose larger set than %d (option --gm-size)\n.", s.gm_size);
	}

	//find the best kappa_data
	set_the_best(ss);

	struct kappa_data *so_far_best = (struct kappa_data*)malloc(sizeof(struct kappa_data));
	kd_init(so_far_best);

	//copy ss->best into so_far_best 
	kd_copy_parameters(ss->best, so_far_best);
	calculate_charges(ss, so_far_best);
	calculate_statistics(ss, so_far_best);

	/* Minimize the result */
	minimize_locally(so_far_best, s.gm_iterations_end);

	/* Tidying up and printing */
	free(bounds);
	kd_copy_parameters(so_far_best, ss->best);
	calculate_charges(ss, ss->best);
	calculate_statistics(ss, ss->best);
	kd_destroy(so_far_best);
	free(so_far_best);
}

/* Generate random population by Latin HyperCube Sampling */
/*void generate_random_population(struct subset* ss, float *bounds, int size) {
	// Get random numbers by Latin Hypercube Sampling
    int dimensions_count = ts.atom_types_count*2+1;
	int points_count = size;
	int seed = rand();
	double* random_lhs = latin_random_new(dimensions_count, points_count, &seed);

	// Redistribute random_lhs[dim_num, point_num] to ss->data 

	//every row contains values of one dimension for all population members
	//start with alpha and beta for all atom types
	for (int i = 0; i < points_count; i++) {
		for (int j = 0; j < ts.atom_types_count; j++) {
			ss->data[i].parameters_alpha[j] = interpolate_to_different_bounds((float)random_lhs[j*points_count + i], bounds[2+j*4], bounds[2+j*4+1]);
			ss->data[i].parameters_beta[j] = interpolate_to_different_bounds((float)random_lhs[(j+1)*points_count + i], bounds[2+j*4+2], bounds[2+j*4+3]);
		}
	}

	//the last row of random_lhs is kappa for all population members
	for (int i = 0; i < points_count; i++) {
		if (s.fixed_kappa < 0)
			ss->data[i].kappa = interpolate_to_different_bounds((float)random_lhs[(dimensions_count-2)*points_count + i], bounds[0], bounds[1]);
		else
			ss->data[i].kappa = s.fixed_kappa;
	}
	free(random_lhs);

}     */

/* Run local minimization on part of population */
int minimize_part_of_gm_set(struct subset* ss, int min_iterations) {
	int quite_good = 0;
	int i = 0;
	//we minimize all with R2>0.2 && R>0
#pragma omp parallel for num_threads(s.gm_threads) shared(ss, quite_good) private(i)
	for (i = 0; i < ss->kappa_data_count; i++) {
		if (ss->data[i].full_stats.R2 > 0.2 && ss->data[i].full_stats.R > 0) {
#pragma omp critical
			{
				quite_good++;
			}
			struct kappa_data* m = (struct kappa_data*) malloc (sizeof(struct kappa_data));
			kd_init(m);
			m->parent_subset = ss;
			kd_copy_parameters(&ss->data[i], m);
			minimize_locally(m, min_iterations);
			kd_copy_parameters(m, &ss->data[i]);
			kd_destroy(m);
			free(m);

		}
	}
	if (s.verbosity >= VERBOSE_KAPPA) {
		printf("Out of %d in population, we minimized %d\n", ss->kappa_data_count, quite_good);
	}
	return quite_good;
}

/* Run local minimization with NEWUOA algorithm */
/*void minimize_locally(struct kappa_data* t, int max_calls) {
	int n = 2*ts.atom_types_count + 1; //number of variables
	int npt = 2*n + 1; //number of interpolation conditions
	double* x = (double*) malloc(n*sizeof(double));
	kappa_data_to_double_array(t, x);
	double rhobeg = 0.2;
	double rhoend = 0.0001;
	int iprint = 0;
	int maxfun = max_calls;
	double* w = (double*) malloc(((npt+13)*(npt+n) + 3*n*(n+3)/2)*sizeof(double));
	//call fortran code NEWUOA for local minimization
	newuoa_(&n, &npt, x, &rhobeg, &rhoend, &iprint, &maxfun, w);
	double_array_to_kappa_data(x, t);
	free(w);
	free(x);
} */

/* Used by NEWUOA algorithm. Evaluates the vector in the local minimization: converts it to kappa_data, computes charges, computes statistics and return the fitness score that should be minimized */
/*extern void calfun_(int n, double*x, double* f) {
	struct kappa_data* t = (struct kappa_data*) malloc (sizeof(struct kappa_data));
	kd_init(t);
	double_array_to_kappa_data(x, t);
	calculate_charges(de_ss, t);
	calculate_statistics_by_sort_mode(t);
	float result = kd_sort_by_return_value(t);
	switch (s.sort_by) {
		case SORT_R:
		case SORT_R2:
		case SORT_RW:
		case SORT_SPEARMAN:
			*f = 1 - (double)(result) + n - n; //+n-n just to get rid of compilaiton warning
			break;
		default:
			*f = (double) (result) + n -n;
	}
	kd_destroy(t);
	free(t);
} */

/* Convert kappa_data into an array of doubles, used in local minimization */
/*void kappa_data_to_double_array(struct kappa_data* t, double* x) {
	x[0] = t->kappa;
	for (int i = 0; i < ts.atom_types_count; i++) {
		x[i*2 + 1] = t->parameters_alpha[i];
		x[i*2 + 2] = t->parameters_beta[i];
	}
} */

/* Convert array of doubles into kappa_data, used in local minimization */
/*void double_array_to_kappa_data(double* x, struct kappa_data* t) {
	t->kappa = (float)x[0];
	for (int i = 0; i < ts.atom_types_count; i++) {
		t->parameters_alpha[i] = (float)x[i*2 + 1];
		t->parameters_beta[i] = (float)x[i*2 + 2];
	}
} */

/* Compute bounds for each parameter of each atom type */
/*void compute_parameters_bounds(float* bounds, int by_atom_type) {
	//returns bounds[kappa_low, kappa_high, alpha_1_low, alpha_1_high, beta_1_low, beta_1_high, alpha_2_low, ...]
	float toH = 0.036749309;
	bounds[0] = 0.0005; //kappa_low
	bounds[1] = 3.5; //kappa_high
	for (int j = 0; j < ts.atom_types_count; j++) {	
		//set bounds to values corresponding with their affinity and ionenergies
		if (by_atom_type == 1) {
			//set bounds for particular atom type
			int atom_number = ts.atom_types[j].Z;
			bounds[2 + j*4] = (((float)(ionenergies[atom_number]) + (float)(affinities[atom_number])*5)/2)*toH; //alpha_low
			bounds[2 + j*4 + 1] = (float) (bounds[2 + j*4] + 0.1); //alpha_high
			bounds[2 + j*4] -= 0.1;
			//bounds[2 + j*4] *= toEV;
			//bounds[2 + j*4 + 1] *= toEV;
			bounds[2 + j*4 + 2] = ((float)(ionenergies[atom_number]) - (float)(affinities[atom_number])*5)/2*toH; //beta_low
			bounds[2 + j*4 + 3]	= (float) (bounds[2 + j*4 + 2] + 0.1); //beta_high
			bounds[2 + j*4 + 2] -= 0.1;
			//bounds[2 + j*4 + 2] *= toEV;
			//bounds[2 + j*4 + 3] *= toEV;
		}
		//set bounds to constant values taken from previous experience
		else if (by_atom_type == 0) {
			bounds[2 + j*4] = 1.8;
			bounds[2 + j*4 + 1] = 3.2;
			bounds[2 + j*4 + 2] = 0;
			bounds[2 + j*4 + 3] = 1.0;
		}
	}
	if (s.verbosity >= VERBOSE_KAPPA) {
		printf("DE Bounds set to:\n");
		for (int i = 0; i < ts.atom_types_count; i++) {
			char buff[10];
			at_format_text(&ts.atom_types[i], buff);
			printf("%s %5.3f - %5.3f, %5.3f - %5.3f\n", buff, bounds[2+i*4], bounds[2+i*4+1], bounds[2+i*4+2], bounds[2+i*4+3]);
		}
	}
} */

/* Get random float within the bounds */
/*float get_random_float(float low, float high) {
	//better random number generator would be nice, but this is sufficient
	float n = low + (float)(rand())/((float)(RAND_MAX/(high-low)));
	return n;
} */

/* Interpolate the number from [0,1] to [low, high] */ 
/*float interpolate_to_different_bounds(float x, float low, float high) {
	return low + x*(high-low);
} */

/* Compute the sum of the vector */
/*int sum(int* vector, int size) {
	int sum = 0;
	for (int i = 0; i < size; i++)
		sum+=vector[i];
	return sum;
} */
