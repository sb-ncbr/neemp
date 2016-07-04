/* Copyright 2013-2016 Tomas Racek (tom@krab1k.net)
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

#ifndef __SETTINGS_H__
#define __SETTINGS_H__

#include <time.h>

#include "config.h"

enum app_mode {

	MODE_PARAMS,
	MODE_CHARGES,
	MODE_INFO,
	MODE_QUALITY,
	MODE_COVER,
	MODE_NOT_SET
};

enum params_calc_method {
	PARAMS_LR_FULL,
	PARAMS_LR_FULL_BRENT,
	PARAMS_DE,
	PARAMS_GM,
	PARAMS_NOT_SET
};

enum sort_mode {

	SORT_R,
	SORT_R2,
	SORT_RW,
	SORT_SPEARMAN,
	SORT_RMSD,
	SORT_RMSD_AVG,
	SORT_D_AVG,
	SORT_D_MAX,
	SORT_NOT_SET
};

enum atom_type_customization {

	AT_CUSTOM_ELEMENT = 0,
	AT_CUSTOM_ELEMENT_BOND = 1,
	AT_CUSTOM_USER = 2
};

enum discarding_mode {

	DISCARD_OFF,
	DISCARD_ITER,
	DISCARD_SIMPLE
};

enum verbosity_levels {

	VERBOSE_MINIMAL = 0,
	VERBOSE_DISCARD = 1,
	VERBOSE_KAPPA = 2
};

struct settings {

	char sdf_file[MAX_PATH_LEN];
	char chg_file[MAX_PATH_LEN];
	char par_file[MAX_PATH_LEN];
	char atb_file[MAX_PATH_LEN];

	char par_out_file[MAX_PATH_LEN];
	char chg_out_file[MAX_PATH_LEN];
	char chg_stats_out_file[MAX_PATH_LEN];

	enum app_mode mode;
	enum params_calc_method params_method;
	enum sort_mode sort_by;
	enum atom_type_customization at_customization;
	enum discarding_mode discard;

	/* Settings regarding PARAMS_LR_FULL* parameters' calculation method */
	float kappa_max;
	float kappa_set;
	float full_scan_precision;

	/* Settings regarding optimization methods as GM, DE, GA */
	int population_size;
	float fixed_kappa; /* Fix kappa to given value */
	int om_threads; /* Number of threads used to paralellize DE */
	int om_iters;
	time_t om_time;
	int polish; /* Use NEWUOA minimization to polish trial or results */

	/* Settings regarding PARAMS_DE method */
	float mutation_constant;
	int dither; /* Set mutation constant to random value from [0.5, 1] each iteration */
	float recombination_constant;

	/* Settings regarding PARAMS_GM optimization method */
	int gm_iterations_beg;
	int gm_iterations_end;

	/* Other settings */
	int random_seed;
	enum verbosity_levels verbosity;

	float tabu_size;

	int limit_iters;
	time_t limit_time;

	int check_charges;
	int list_omitted_molecules;

	int max_threads;

	int extra_precise;
};

void s_init(void);

void parse_options(int argc, char **argv);
void check_settings(void);
void print_settings(void);

char *get_atom_types_by_string(enum atom_type_customization atc);

#endif /* __SETTINGS_H__ */
