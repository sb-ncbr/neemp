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

struct settings {

	char sdf_filename[MAX_PATH_LEN];
	char chg_filename[MAX_PATH_LEN];
	char par_filename[MAX_PATH_LEN];
};

void parse_options(int argc, char **argv);
void check_settings(void);

#endif /* __SETTINGS_H__ */
