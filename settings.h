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

#endif /* __SETTINGS_H__ */
