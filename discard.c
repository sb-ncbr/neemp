/*
 * NEEMP - discard.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#include <stdio.h>
#include <stdlib.h>

#include "discard.h"
#include "kappa.h"
#include "settings.h"
#include "subset.h"
#include "structures.h"

extern struct settings s;
extern struct training_set ts;

/* Perform iterative discarding by traversing the state space */
struct subset *discard_iterative(const struct subset * const initial) {

	struct subset *current, *old;

	/* We guarantee that we won't free subset pointed by intial */
	old = (struct subset *) initial;
	current = old;

	for(int i = 0; i < 3; i++) {
		current = (struct subset *) calloc(1, sizeof(struct subset));
		current->parent = old;

		/* Flip just one molecule from the parent */
		b_init(&current->molecules, ts.molecules_count);
		b_set_as(&current->molecules, &old->molecules);
		int bno = rand() % ts.molecules_count;
		b_flip(&current->molecules, bno);

		if(s.verbosity >= VERBOSE_DISCARD) {
			fprintf(stdout, "\nIteration no. %d. Molecule discarded: %s Molecules used: %d\n",\
				i, ts.molecules[bno].name, b_count_bits(&current->molecules));
		}

		/* The actual calculations */
		find_the_best_parameters_for_subset(current);

		if(s.verbosity >= VERBOSE_DISCARD) {
			printf("=> ");
			kd_print_stats(current->best);
			printf("\n");
		}

		/* Check if the new subset is better than the old one */
		if(kd_sort_by_is_better(current->best, old->best)) {

			/* Do not free initial subset! */
			if(old != initial) {
				ss_destroy(old);
				free(old);
			}

			old = current;
		}
		else {
			ss_destroy(current);
			free(current);

			/* If this is the last iteration, we want to return the best result so far */
			current = old;
		}
	}

	return current;
}

/* Perform simple discard in the way that EMP does */
struct subset *discard_simple(const struct subset * const initial) {

	struct subset *best = (struct subset *) initial;
	struct subset *current;

	for(int i = 0; i < ts.molecules_count && i < 3; i++) {

		current = (struct subset *) calloc(1, sizeof(struct subset));
		current->parent = best;

		/* Flip i-th molecule from the parent */
		b_init(&current->molecules, ts.molecules_count);
		b_set_as(&current->molecules, &best->molecules);
		b_flip(&current->molecules, i);

		if(s.verbosity >= VERBOSE_DISCARD) {
			fprintf(stdout, "\nIteration no. %d. Molecule discarded: %s Molecules used: %d\n",\
				i, ts.molecules[i].name, b_count_bits(&current->molecules));
		}

		/* The actual calculations */
		find_the_best_parameters_for_subset(current);

		if(s.verbosity >= VERBOSE_DISCARD) {
			printf("=> ");
			kd_print_stats(current->best);
			printf("\n");
		}

		if(kd_sort_by_is_better(current->best, best->best)) {

			/* Do not free initial subset! */
			if(best != initial) {
				ss_destroy(best);
				free(best);
			}

			best = current;
		}
		else {
			ss_destroy(current);
			free(current);
		}
	}

	return best;
}
