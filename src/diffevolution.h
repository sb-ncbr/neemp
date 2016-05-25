/*
 * NEEMP - diffevolution.h
 *
 * by Jana Pazurikova (pazurikova@ics.muni.cz)
 * 2016
 *
 * */

#ifndef __DIFFEVOLUTION_H__
#define __DIFFEVOLUTION_H__

#include "subset.h"

void run_diff_evolution(struct subset * const ss);
void generate_random_population(struct subset* ss, float *bounds, int size);
int minimize_part_of_population(struct subset* ss, int* good_indices);
int evolve_kappa(struct kappa_data* trial, struct kappa_data* x, struct kappa_data* a, struct kappa_data *b, float *bounds, float mutation_constant, float recombination_constant);
int compare_and_set(struct kappa_data* trial, struct kappa_data* so_far_best);
void compute_parameters_bounds(float* bounds, int by_atom_type);
float get_random_float(float low, float high);
float interpolate_to_different_bounds(float x, float low, float high);
int sum(int* vector, int size);
void kd_copy_parameters(struct kappa_data* from, struct kappa_data* to);
void kd_copy_statistics(struct kappa_data* from, struct kappa_data* to);
void minimize_locally(struct kappa_data* trial, int max_calls);
extern void newuoa_(int* n, int* npt, double* x, double* rhobeg, double* rhoend, int* iprint, int* maxfun, double* w);
extern void calfun_(int n, double*x, double* f);
void kappa_data_to_double_array(struct kappa_data* trial, double* x);
void double_array_to_kappa_data(double* x, struct kappa_data* trial);
int is_quite_good(struct kappa_data* t);

struct subset * de_ss;
#endif /* __DIFFEVOLUTION_H__ */
