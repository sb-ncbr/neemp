/*
 * NEEMP - settings.h
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013, 2014
 *
 * */

#ifndef __SETTINGS_H__
#define __SETTINGS_H__

#include <time.h>

#include "config.h"

enum app_mode {

	MODE_PARAMS,
	MODE_CHARGES,
	MODE_INFO,
	MODE_CROSS,
	MODE_COVER,
	MODE_NOT_SET
};

enum params_calc_method {
	PARAMS_LR_FULL,
	PARAMS_LR_FULL_BRENT,
	PARAMS_DE,
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
	SORT_D_MAX
};

enum atom_type_customization {

	AT_CUSTOM_ELEMENT = 0,
	AT_CUSTOM_ELEMENT_BOND = 1,
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
	char wgh_file[MAX_PATH_LEN];

	char par_out_file[MAX_PATH_LEN];
	char chg_out_file[MAX_PATH_LEN];
	char chg_stats_out_file[MAX_PATH_LEN];

	enum app_mode mode;
	enum params_calc_method params_method;
	enum sort_mode sort_by;
	enum atom_type_customization at_customization;
	enum discarding_mode discard;

	//settings regarding PARAMS_LR_FULL* parameters' calculation method
	float kappa_max;
	float kappa_set;
	float full_scan_precision;

	//settings regarding PARAMS_DE optimization method
	int population_size;
	float mutation_constant;
	int dither; //set mutation constant to random value from [0.5, 1] each iteration
	float recombination_constant;
	float fixed_kappa; //fix kappa to given value
	int de_threads; //number of threads used to paralellize DE
	int limit_de_iters;
	time_t limit_de_time;
	int polish; //use NEWUOA minimization to polish trial or results

	//other settings
	int random_seed;
	enum verbosity_levels verbosity;

	float tabu_size;

	int limit_iters;
	time_t limit_time;

	int check_charges;
	int list_omitted_molecules;

	int max_threads;

	float rw;
};

void s_init(void);

void parse_options(int argc, char **argv);
void check_settings(void);
void print_settings(void);

char *get_atom_types_by_string(enum atom_type_customization atc);

#endif /* __SETTINGS_H__ */
