/*
 * NEEMP - guidedmin.h
 *
 * by Jana Pazurikova (pazurikova@ics.muni.cz)
 * 2016
 *
 * */

#ifndef __GUIDEDMIN_H__
#define __GUIDEDMIN_H__

#include "subset.h"
#include "diffevolution.h"

void run_guided_min(struct subset * const ss);
//void generate_random_population(struct subset* ss, float *bounds, int size);
int minimize_part_of_gm_set(struct subset* ss, int min_iterations);
//void compute_parameters_bounds(float* bounds, int by_atom_type);
//float get_random_float(float low, float high);
//float interpolate_to_different_bounds(float x, float low, float high);
//int sum(int* vector, int size);
//void minimize_locally(struct kappa_data* trial, int max_calls);
//extern void newuoa_(int* n, int* npt, double* x, double* rhobeg, double* rhoend, int* iprint, int* maxfun, double* w);
//extern void calfun_(int n, double*x, double* f);
//void kappa_data_to_double_array(struct kappa_data* trial, double* x);
//void double_array_to_kappa_data(double* x, struct kappa_data* trial);

struct subset * de_ss;
#endif /* __GUIDEDMIN_H__ */
