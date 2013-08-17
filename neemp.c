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

	s_init();

	parse_options(argc, argv);
	check_settings();

	if(s.mode == MODE_PARAMS) {
		load_molecules();
		load_charges();
		preprocess_molecules();

	} else if(s.mode == MODE_CHARGES) {
		load_molecules();
		load_parameters();

	} else if(s.mode == MODE_INFO) {
		load_molecules();
	}

	ts_destroy();

	return RETURN_OK;
}
