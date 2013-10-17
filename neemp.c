/*
 * NEEMP - neemp.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#include <mkl.h>
#include <string.h>

#include "eem.h"
#include "kappa.h"
#include "neemp.h"
#include "io.h"
#include "settings.h"
#include "subset.h"
#include "statistics.h"
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
		b_init(&full.molecules, ts.molecules_count);
		b_set_all(&full.molecules);
		find_the_best_parameters_for_subset(&full);

		if(s.chg_stats_out_file[0] != '\0')
			output_charges_stats(&full);

		if(s.par_out_file[0] != '\0')
			output_parameters(&full);

		print_results(&full);

		ss_destroy(&full);

	} else if(s.mode == MODE_CHARGES) {
		load_molecules();

		preprocess_molecules();

		struct subset full;
		b_init(&full.molecules, ts.molecules_count);
		b_set_all(&full.molecules);

		full.kappa_data_count = 1;
		full.data = (struct kappa_data *) calloc(1, sizeof(struct kappa_data));
		full.best = &full.data[0];
		kd_init(full.best);

		load_parameters(full.best);

		calculate_charges(&full, full.best);
		calculate_statistics(&full, full.best);

		print_results(&full);

		output_charges(&full);

		ss_destroy(&full);

	} else if(s.mode == MODE_INFO) {
		load_molecules();
		preprocess_molecules();
		ts_info();
	}

	ts_destroy();

	mkl_free_buffers();

	return RETURN_OK;
}
