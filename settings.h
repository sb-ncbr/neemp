/*
 * NEEMP - settings.h
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#ifndef __SETTINGS_H__
#define __SETTINGS_H__

#include "config.h"

enum app_mode {
	MODE_PARAMS,
	MODE_CHARGES,
	MODE_INFO,
	MODE_NOT_SET
};

struct settings {

	char sdf_filename[MAX_PATH_LEN];
	char chg_filename[MAX_PATH_LEN];
	char par_filename[MAX_PATH_LEN];
	char chgout_filename[MAX_PATH_LEN];

	enum app_mode mode;
	int full_scan_only;

	float kappa_max;
	float kappa_set;
	float full_scan_precision;
};

void s_init(void);

void parse_options(int argc, char **argv);
void check_settings(void);

#endif /* __SETTINGS_H__ */
