/* Copyright 2013-2016 Jana Pazurikova (pazurikova@ics.muni.cz)
 *
 * This file is part of NEEMP.
 *
 * NEEMP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * NEEMP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with NEEMP. If not, see <http://www.gnu.org/licenses/>.
 */

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
void minimize_locally(struct kappa_data* trial, int max_calls);
extern void newuoa_(int* n, int* npt, double* x, double* rhobeg, double* rhoend, int* iprint, int* maxfun, double* w);
extern void calfun_(int n, double*x, double* f);
void kappa_data_to_double_array(struct kappa_data* trial, double* x);
void double_array_to_kappa_data(double* x, struct kappa_data* trial);
int is_quite_good(const struct kappa_data * const t);

struct subset * de_ss;
#endif /* __DIFFEVOLUTION_H__ */
