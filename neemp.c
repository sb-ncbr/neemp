/*
 * NEEMP - neemp.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#include <mkl.h>
#include <string.h>

#include "kappa.h"
#include "neemp.h"
#include "io.h"
#include "settings.h"
#include "subset.h"
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

		ts_info();

		struct subset full;
		find_the_best_parameters_for_subset(&full);

		if(s.chgout_filename[0] != '\0')
			output_charges_stats(&full);

		ss_destroy(&full);

	} else if(s.mode == MODE_CHARGES) {
		load_molecules();
		load_parameters();

	} else if(s.mode == MODE_INFO) {
		load_molecules();
		preprocess_molecules();
		ts_info();
	}

	ts_destroy();

	mkl_free_buffers();

	return RETURN_OK;
}
