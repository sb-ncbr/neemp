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
#include "diffevolution.h"

extern const struct training_set ts;
extern const struct settings s;
extern const float ionenergies[];
extern const float affinities[];


/* Run differential evolution algorithm to find the best set of parameters for calculation of partial charges. */ 
void run_diff_evolution(struct subset * const ss) {

	/* Create a set of random points in vector space of kappa_data */            
	if (s.verbosity >= VERBOSE_KAPPA)
		printf("DE Generating population of size %d\n", s.population_size);

	/* Set bounds for each parameters in kappa_data */
	float *bounds = (float*) malloc((ts.atom_types_count*2+1)*2*sizeof(float));
	//compute bounds, 0 means set them to fixed numbers taken from Tomas's full scan, 1 means try to find them with broad search
	compute_parameters_bounds(ss,bounds, 0);
	//we must create "regular" population only after computing bounds, as sometimes we can compute bounds with usage of preliminary population
	fill_ss(ss, s.population_size); 
	generate_random_population(ss, bounds, s.population_size);
	de_ss = ss;
	if (s.polish > 2) {
		if (s.verbosity >= VERBOSE_KAPPA)
			printf("DE minimizing part of population\n");
		minimize_part_of_population(ss, (int)floor(s.population_size*0.25));
	}
	/* Evaluate the fitness function for all points and assign the best */
	if (s.verbosity >= VERBOSE_KAPPA)
		printf("DE Calculating charges and evaluating fitness function for whole population\n");
	int i = 0;
#pragma omp parallel for num_threads(s.de_threads) default(shared) private(i)
	for (i = 0; i < ss->kappa_data_count; i++) {
		calculate_charges(ss, &ss->data[i]);
		calculate_statistics_by_sort_mode(&ss->data[i]);
	}
	//TODO extract to separate method, also used in kappa.c:find_the_best_parameters
	ss->best = &ss->data[0];
	for (i = 0; i < ss->kappa_data_count -1; i++)
		if (kd_sort_by_is_better(&ss->data[i], ss->best))
			ss->best = &ss->data[i];

	/* Run the optimization for max iterations */
	//TODO include iters_max for DE in limits or use one already there (that is used for discard)
	//TODO can we tell we have converged? if yes, include to while condition
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
	kd_print_results(so_far_best);
	minimize_locally(so_far_best, 500);
	calculate_charges(ss, so_far_best);
	calculate_statistics(ss, so_far_best);
	kd_print_results(so_far_best);
	float mutation_constant = s.mutation_constant;
	int iters_with_evolution=0;
	int condition=1;

#pragma omp parallel num_threads(s.de_threads) default(shared)
	{
#pragma omp master
		while (condition) {
			condition = ((iter < s.limit_de_iters) || ((iter < 2*s.limit_de_iters) && (!is_good_enough(so_far_best))));
			{
#pragma omp master
				{
					iter++;
					if (s.verbosity >= VERBOSE_KAPPA) {
						if (iter % 100 == 0)
							printf("\nDE iter %d evolve thread %d\n", iter, omp_get_thread_num());
						else
							printf(".");
					}
					/* Select randomly two points from population */
					//TODO replace with some real random number generator
					int rand1 = (int)(floor(get_random_float(0, (float)ss->kappa_data_count -1 )));
					int rand2 = (int)(floor(get_random_float(0, (float)ss->kappa_data_count -1 )));

					struct kappa_data* a = &(ss->data[rand1]);
					struct kappa_data* b = &(ss->data[rand2]);

					if (s.dither)
						mutation_constant = get_random_float(0.5, 1);
					/* Recombine parts of best, a and b to obtain new trial structure */
#pragma omp critical
					evolve_kappa(trial, so_far_best, a, b, bounds, mutation_constant, s.recombination_constant);
					iters_with_evolution++;
					/* Evaluate the new trial structure */
					calculate_charges(ss, trial);
					calculate_statistics_by_sort_mode(trial);
					/* If the new structure is better than what we have before, reassign */
#pragma omp critical
					{
						if (compare_and_set(trial, so_far_best)) {
							calculate_charges(ss, so_far_best);
							calculate_statistics(ss, so_far_best);
							if (s.verbosity >= VERBOSE_KAPPA) {
								printf("\n");
								kd_print_results(so_far_best);
							}
						}
					}
				}
				/* All other threads do this */ 
				if (s.polish > 1 && (is_quite_good(trial) || (s.de_threads == 0 || omp_get_thread_num() != 0)))
				{
					if (s.verbosity > VERBOSE_KAPPA)
						printf("\nDE min thread %d\n", omp_get_thread_num());
					struct kappa_data *min_trial = (struct kappa_data*)malloc(sizeof(struct kappa_data));
					kd_init(min_trial);
					min_trial->parent_subset = ss;
					kd_copy_parameters(trial, min_trial);
					minimize_locally(min_trial, 100);
					calculate_charges(de_ss, min_trial);
					calculate_statistics(de_ss, min_trial);
#pragma omp critical
					{
						if (kd_sort_by_is_better(min_trial, trial) && compare_and_set(min_trial, so_far_best)) {
							calculate_charges(ss, so_far_best);
							calculate_statistics(ss, so_far_best);
							if(s.verbosity >= VERBOSE_KAPPA) {
								printf("\n");
								kd_print_results(so_far_best);
							}
						}
					}
					kd_destroy(min_trial);
					free(min_trial);
				}
			}
		}
	}
	if (s.polish > 0)
		minimize_locally(so_far_best, 1000);
	kd_destroy(trial);
	free(trial);
	free(bounds);
	kd_copy_parameters(so_far_best, ss->best);
	calculate_charges(ss, ss->best);
	calculate_statistics(ss, ss->best);
}

/* Generate random population by Latin HyperCube Sampling */
void generate_random_population(struct subset* ss, float *bounds, int size) {
	/* Get random numbers by Latin Hypercube Sampling */
	int dimensions_count = ts.atom_types_count*2+1;
	int points_count = size;
	int seed = rand(); //get_seed();
	double* random_lhs = latin_random_new(dimensions_count, points_count, &seed);

	/* Redistribute random_lhs[dim_num, point_num] to ss->data */

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

}

void minimize_part_of_population(struct subset* ss, int count) {
	struct kappa_data* m = (struct kappa_data*) malloc (sizeof(struct kappa_data));
	kd_init(m);
	m->parent_subset = ss;
	for (int i = 0; i < count; i++) {
		int r = (int) (floor(get_random_float(0, (float)ss->kappa_data_count-1)));
		kd_copy_parameters(&ss->data[r], m);
		minimize_locally(m, 500);
		kd_copy_parameters(m, &ss->data[r]);

	}
	kd_destroy(m);
	free(m);
}

/* Evolve kappa_data, i.e. create a new trial structure */
int evolve_kappa(struct kappa_data* trial, struct kappa_data* x, struct kappa_data* a, struct kappa_data *b, float *bounds, float mutation_constant, float recombination_constant) {
	int changed = 0;
	kd_copy_parameters(x, trial);
	for (int i = 0; i < ts.atom_types_count; i++)
		if (get_random_float(0,1) < recombination_constant)	{
			changed++;
			trial->parameters_alpha[i] += mutation_constant*(a->parameters_alpha[i] - b->parameters_alpha[i]);
			if (bounds[2 + i*4] > trial->parameters_alpha[i] || bounds[2 + i*4 + 1] < trial->parameters_alpha[i]) {
				trial->parameters_alpha[i] = x->parameters_alpha[i];
				changed--;
			}
		}
	for (int i = 0; i < ts.atom_types_count; i++)
		if (get_random_float(0,1) < recombination_constant)	{
			changed++;
			trial->parameters_beta[i] += mutation_constant*(a->parameters_beta[i] - b->parameters_beta[i]);
			if (bounds[2 + i*4 + 2] > trial->parameters_beta[i] || bounds[2 + i*4 + 3] < trial->parameters_beta[i]) {
				trial->parameters_beta[i] = x->parameters_beta[i];
				changed--;
			}
		}

	//always change kappa
	if (s.fixed_kappa < 0)	{
		trial->kappa += mutation_constant*(a->kappa - b->kappa)/(bounds[1]-bounds[0]);
		changed++;
		if (bounds[0] > trial->kappa || bounds[1] < trial->kappa) {
			trial->kappa -= mutation_constant*(a->kappa - b->kappa)/(bounds[1]-bounds[0]);
		}
	}
	else 
		trial->kappa = s.fixed_kappa;
	return changed;
}

/* Compare trial structure with so far best structure and save the better */
int compare_and_set(struct kappa_data* trial, struct kappa_data* so_far_best) {
	/* If we consider only total statistics */
	if (s.evolve_by_element == 0) {
		if (kd_sort_by_is_better(trial, so_far_best)) {
			kd_copy_parameters(trial, so_far_best);
			return 1;
		}
		else
			return 0;
	}

	/* We consider also per atom statistics */
	float compare_full_stats = kd_sort_by_return_value(trial) - kd_sort_by_return_value(so_far_best);
	int trial_is_better = compare_full_stats > 0;
	compare_full_stats = fabsf(compare_full_stats);
	float compare_kappas = fabsf(trial->kappa - so_far_best->kappa);
	float kappas_threshold = 0.05;
	float full_stats_threshold = 0.01;
	float per_atom_threshold = 0.1;
	//if full_stats are much different, we change or not change according to them
	if (compare_full_stats >= full_stats_threshold) {
		if (trial_is_better) {
			kd_copy_parameters(trial, so_far_best);
			return 1;
		}
		else 
			return 0;
	}

	//if kappas are too different, we behave according to full_stats
	if (compare_kappas >= kappas_threshold) {
		if (trial_is_better) {
			kd_copy_parameters(trial, so_far_best);
			return 1;
		}
		else
			return 0;
	}

	//if we get here, we might consider change to a bit worse or keep a bit worse if per atom stats show great differences
	int* better_per_atom = (int*) malloc(ts.atom_types_count*sizeof(int));
	kd_sort_by_is_much_better_per_atom(better_per_atom, trial, so_far_best, per_atom_threshold);
	//if kappa is close, trial is a bit worse in full stats, but much better in one (or more) of the elements, than allow a bit of worsening and set so_far_best to trial 
	if (!trial_is_better && sum(better_per_atom, ts.atom_types_count)) {
		if (s.verbosity >= VERBOSE_KAPPA)
			printf("DE allowing a bit worse\n");
		kd_copy_parameters(trial, so_far_best);
		return 1;
	}

	//if trial is a bit better in full stats, but much worse in one of the elements, that prevent a bettering and keep so_far_best the same
	if (trial_is_better && sum(better_per_atom, ts.atom_types_count) < 0) {
		if (s.verbosity >= VERBOSE_KAPPA)
			printf("DE preventing a bit better\n");
		return 0;
	}
	return 0;
}

void minimize_locally(struct kappa_data* t, int max_calls) {
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
}

extern void calfun_(int n, double*x, double* f) {
	struct kappa_data* t = (struct kappa_data*) malloc (sizeof(struct kappa_data));
	kd_init(t);
	double_array_to_kappa_data(x, t);
	calculate_charges(de_ss, t);
	calculate_statistics_by_sort_mode(t);
	float result = kd_sort_by_return_value(t);
	*f = 1 - (double)(result) + n - n; //+n-n just to get rid of compilaiton warning
	kd_destroy(t);
	free(t);
}

void kappa_data_to_double_array(struct kappa_data* t, double* x) {
	x[0] = t->kappa;
	for (int i = 0; i < ts.atom_types_count; i++) {
		x[i*2 + 1] = t->parameters_alpha[i];
		x[i*2 + 2] = t->parameters_beta[i];
	}
}

void double_array_to_kappa_data(double* x, struct kappa_data* t) {
	t->kappa = (float)x[0];
	for (int i = 0; i < ts.atom_types_count; i++) {
		t->parameters_alpha[i] = (float)x[i*2 + 1];
		t->parameters_beta[i] = (float)x[i*2 + 2];
	}
}

int is_good_enough(struct kappa_data* t) {
	if (t->full_stats.R2 < 0.85)
		return 0;
	for (int i = 0; i < ts.atom_types_count; i++) {
		if (t->per_at_stats[i].R2 < 0.6)
			return 0; 
	}
	return 1;
}

int is_quite_good(struct kappa_data* t) {
	if (t->full_stats.R2 < 0.6)
		return 0;
	for (int i = 0; i < ts.atom_types_count; i++) {
		if (t->per_at_stats[i].R2 < 0.5)
			return 0; 
	}
	return 1;
}

/* Compute bounds for each parameter of each atom type */
void compute_parameters_bounds(struct subset *ss, float* bounds, int by_atom_type) {          
	//returns bounds[kappa_low, kappa_high, alpha_1_low, alpha_1_high, beta_1_low, beta_1_high, alpha_2_low, ...]
	float toH = 0.036749309;
	bounds[0] = 0.0005; //kappa_low
	bounds[1] = 3.5; //kappa_high
	for (int j = 0; j < ts.atom_types_count; j++) {	
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
		else if (by_atom_type == 0) {
			bounds[2 + j*4] = 2;
			bounds[2 + j*4 + 1] = 3;
			bounds[2 + j*4 + 2] = 0;//0;
			bounds[2 + j*4 + 3] = 0.8;
		}
		else if (by_atom_type == 2) {
			bounds[2 + j*4] = 0;
			bounds[2 + j*4 + 1] = 20;
			bounds[2 + j*4 + 2] = 0;
			bounds[2 + j*4 + 3] = 5;
		}
	}

	/* Restrict bounds for each atom type to reasonable values by one iteration of DE to find out what values gives the best results */
	if (by_atom_type == 2) {
		int restrict_bounds_population_size = 40;
		fill_ss(ss, restrict_bounds_population_size);
		generate_random_population(ss, bounds, restrict_bounds_population_size);
		struct kappa_data** best_per_atom = (struct kappa_data**) malloc((ts.atom_types_count)*(sizeof(struct kappa_data*)));
		for (int i = 0; i < restrict_bounds_population_size; i++) {
			calculate_charges(ss, &ss->data[i]);
			calculate_statistics(ss, &ss->data[i]);
		}

		for (int i = 0; i < ts.atom_types_count; i++) 
			best_per_atom[i] = &ss->data[0];

		for (int i = 0; i < ts.atom_types_count; i++) {
			for (int j = 0; j< restrict_bounds_population_size; j++) {
				if (kd_sort_by_is_better_per_atom(&ss->data[j], best_per_atom[i], i))
					best_per_atom[i] = &ss->data[j];
			}
		}

		for (int i = 0; i < ts.atom_types_count; i++) {
			char buff[10];
			at_format_text(&ts.atom_types[i], buff);
			printf("Best for atom %s with %f\n", buff, best_per_atom[i]->per_at_stats[i].R);
			print_parameters(best_per_atom[i]);
			if (best_per_atom[i]->per_at_stats[i].R > 0.2)
			{
				bounds[2 + i*4] = (float)0.8*(best_per_atom[i]->parameters_alpha[i]);
				bounds[2 + i*4 + 1] = (float)1.2*(best_per_atom[i]->parameters_alpha[i]);
				bounds[2 + i*4 + 2] = (float)0.8*(best_per_atom[i]->parameters_beta[i]);
				bounds[2 + i*4 + 3] = (float)1.2*(best_per_atom[i]->parameters_beta[i]);
			}
		}

		for (int i = 0; i < restrict_bounds_population_size; i++) {
			kd_destroy(&ss->data[i]);
		}
		free(ss->data);


	}

	if (s.verbosity >= VERBOSE_KAPPA) {
		printf("DE Bounds set to:\n");
		for (int i = 0; i < ts.atom_types_count; i++) {
			char buff[10];
			at_format_text(&ts.atom_types[i], buff);
			printf("%s %5.3f - %5.3f, %5.3f - %5.3f\n", buff, bounds[2+i*4], bounds[2+i*4+1], bounds[2+i*4+2], bounds[2+i*4+3]);
		}
	}
}

/* Get random float within the bounds */
float get_random_float(float low, float high) {
	float n = low + (float)(rand())/((float)(RAND_MAX/(high-low)));
	return n;
}

/* Interpolate the number from [0,1] to [low, high] */ 
float interpolate_to_different_bounds(float x, float low, float high) {
	return low + x*(high-low);
}

/* Compute the sum of the vector */
int sum(int* vector, int size) {
	int sum = 0;
	for (int i = 0; i < size; i++)
		sum+=vector[i];
	return sum;
}

/* Copy data from one kappa_data to another */
void kd_copy_parameters(struct kappa_data* from, struct kappa_data* to) {
	to->kappa = from->kappa;
	for (int i = 0; i < ts.atom_types_count; i++) {
		to->parameters_alpha[i] = from->parameters_alpha[i];
		to->parameters_beta[i] = from->parameters_beta[i];
	}
}
