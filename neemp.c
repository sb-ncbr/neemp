/*
 * NEEMP - neemp.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#include <assert.h>
#include <signal.h>
#include <string.h>

#ifdef USE_MKL
#include <mkl.h>
#endif /* USE_MKL */

#include "discard.h"
#include "eem.h"
#include "kappa.h"
#include "neemp.h"
#include "io.h"
#include "limits.h"
#include "settings.h"
#include "subset.h"
#include "statistics.h"
#include "structures.h"

struct training_set ts;
struct settings s;
struct limit limits;

int termination_flag = 0;

static void sig_handler(int sig __attribute__ ((unused)));

/* Set termination flag on signal received */
static void sig_handler(int sig __attribute__ ((unused))) {

	termination_flag = 1;
}

/* That's the main thing */
int main(int argc, char **argv) {

	time_t start, end;

	start = time(NULL);

	s_init();
	parse_options(argc, argv);
	check_settings();

	printf("%s (%s) started\n", APP_NAME, APP_VERSION);

	l_init(&limits, s.limit_iters, s.limit_time);

	load_molecules();

	/* Interrupt discarding if one of these signals is received */
	signal(SIGINT, sig_handler);
	signal(SIGTERM, sig_handler);

	switch(s.mode) {
		case MODE_PARAMS: {
			load_charges();
			preprocess_molecules();
			ts_info();

			struct subset full;
			b_init(&full.molecules, ts.molecules_count);
			b_set_all(&full.molecules);
			find_the_best_parameters_for_subset(&full);

			printf("\nResults for the full set:\n\n");
			print_results(&full);
			struct subset *result = NULL;

			switch(s.discard) {
				case DISCARD_OFF:
					result = &full;
					break;
				case DISCARD_ITER: {
					result = discard_iterative(&full);
					break;
				}
				case DISCARD_SIMPLE: {
					result = discard_simple(&full);
					break;
				}
			}

			if(s.chg_stats_out_file[0] != '\0')
				output_charges_stats(result);

			if(s.par_out_file[0] != '\0')
				output_parameters(result);

			if(result != &full) {
				printf("\nFinal results after discarding:\n\n");
				print_results(result);

				/* Clean up the discarding result */
				ss_destroy(result);
				free(result);
			}
			ss_destroy(&full);
			break;
		}
		case MODE_CHARGES: {
			struct subset full;
			preprocess_molecules();
			b_init(&full.molecules, ts.molecules_count);
			b_set_all(&full.molecules);

			/* We use only one particular kappa which is read from .par file */
			full.kappa_data_count = 1;
			full.data = (struct kappa_data *) calloc(1, sizeof(struct kappa_data));
			full.best = &full.data[0];
			kd_init(full.best);

			load_parameters(full.best);

			calculate_charges(&full, full.best);
			output_charges(&full);

			ss_destroy(&full);
			break;
		}
		case MODE_INFO:
			preprocess_molecules();
			ts_info();
			break;
		case MODE_NOT_SET:
			/* Something bad happened. */
			assert(0);
	}

	ts_destroy();

	#ifdef USE_MKL
	mkl_free_buffers();
	#endif /* USE_MKL */

	end = time(NULL);
	time_t diff = end - start;
	time_t hours = diff / 3600;
	diff %= 3600;
	time_t minutes = diff / 60;
	diff %= 60;

	printf("Execution took %lu hour(s) %lu minute(s) %lu second(s)\n",\
		 (unsigned long) hours, (unsigned long) minutes, (unsigned long) diff);
	printf("%s (%s) ended.\n", APP_NAME, APP_VERSION);

	return RETURN_OK;
}
