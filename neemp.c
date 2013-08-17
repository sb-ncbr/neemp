/*
 * NEEMP - neemp.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#include <string.h>

#include "neemp.h"
#include "parser.h"
#include "settings.h"
#include "structures.h"

struct training_set ts;
struct settings s;

int main(int argc, char **argv) {

	parse_options(argc, argv);
	check_settings();

	load_molecules();
	load_charges();

	preprocess_molecules();

	ts_destroy();

	return RETURN_OK;
}
