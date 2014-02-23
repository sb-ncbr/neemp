/*
 * NEEMP - settings.h
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
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
	MODE_NOT_SET
};

enum sort_mode {

	SORT_R,
	SORT_RMSD,
	SORT_MSE,
	SORT_D_AVG,
	SORT_D_MAX
};

enum atom_type_customization {

	AT_CUSTOM_ELEMENT,
	AT_CUSTOM_ELEMENT_BOND,
	AT_CUSTOM_PARTNER,
	AT_CUSTOM_VALENCE
};

enum discarding_mode {

	DISCARD_OFF,
	DISCARD_ITER,
	DISCARD_SIMPLE
};

enum verbosity_levels {

	VERBOSE_DISCARD = 1,
	VERBOSE_KAPPA = 2
};

struct settings {

	char sdf_file[MAX_PATH_LEN];
	char chg_file[MAX_PATH_LEN];
	char par_file[MAX_PATH_LEN];

	char par_out_file[MAX_PATH_LEN];
	char chg_out_file[MAX_PATH_LEN];
	char chg_stats_out_file[MAX_PATH_LEN];

	enum app_mode mode;
	enum sort_mode sort_by;
	enum atom_type_customization at_customization;
	enum discarding_mode discard;
	int full_scan_only;

	float kappa_max;
	float kappa_set;
	float full_scan_precision;

	enum verbosity_levels verbosity;

	float tabu_size;

	int limit_iters;
	time_t limit_time;
};

void s_init(void);

void parse_options(int argc, char **argv);
void check_settings(void);
void print_settings(void);

#endif /* __SETTINGS_H__ */
