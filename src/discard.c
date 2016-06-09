/*
 * NEEMP - discard.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013, 2014
 *
 * */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "discard.h"
#include "kappa.h"
#include "limits.h"
#include "settings.h"
#include "subset.h"
#include "structures.h"
#include "tabu.h"

extern struct settings s;
extern struct training_set ts;
extern struct limit limits;
extern int termination_flag;

/* Perform iterative discarding by traversing the state space */
struct subset *discard_iterative(const struct subset * const initial) {

	assert(initial != NULL);

	struct subset *current, *old;

	/* We guarantee that we won't free subset pointed by intial */
	old = (struct subset *) initial;
	current = old;

	/* Initialize tabu list */
	struct tabu ban_list;
	t_init(&ban_list, ts.molecules_count * s.tabu_size);

	for(int i = 0; !l_check(&limits) && !termination_flag; i++) {
		/* Flip just one molecule from the parent */
		current = (struct subset *) calloc(1, sizeof(struct subset));
		ss_init(current, old);

		/* Generate random molecule to flip, check if it's allowed */
		int mol_idx;
		do {
			mol_idx = rand() % ts.molecules_count;
		} while(t_is_banned(&ban_list, mol_idx));
		t_update(&ban_list, mol_idx);

		b_flip(&current->molecules, mol_idx);

		if(s.verbosity >= VERBOSE_DISCARD) {
			fprintf(stdout, "\nIteration no. %d. Molecule: %s Molecules used: %d\n",\
				i + 1, ts.molecules[mol_idx].name, b_count_bits(&current->molecules));
		}

		/* The actual calculations */
		find_the_best_parameters_for_subset(current);

		if(s.verbosity >= VERBOSE_DISCARD) {
			printf(" %c ", kd_sort_by_is_better(current->best, old->best) ? '+' : '-');
			kd_print_stats(current->best);
			print_results(current);
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

		/* Update number of iterations */
		limits.iters_current += 1;
	}

	t_destroy(&ban_list);

	return current;
}

/* Perform simple discard in the way that EMP does */
struct subset *discard_simple(const struct subset * const initial) {

	assert(initial != NULL);

	struct subset *best = (struct subset *) initial;
	struct subset *current;

	for(int i = 0; i < ts.molecules_count && !termination_flag && !l_check(&limits); i++) {

		current = (struct subset *) calloc(1, sizeof(struct subset));
		ss_init(current, best);

		/* Flip i-th molecule from the parent */
		b_flip(&current->molecules, i);

		if(s.verbosity >= VERBOSE_DISCARD) {
			fprintf(stdout, "\nIteration no. %d. Molecule: %s Molecules used: %d\n",\
				i + 1, ts.molecules[i].name, b_count_bits(&current->molecules));
		}

		/* The actual calculations */
		find_the_best_parameters_for_subset(current);

		if(s.verbosity >= VERBOSE_DISCARD) {
			printf(" %c ", kd_sort_by_is_better(current->best, best->best) ? '+' : '-');
			kd_print_stats(current->best);
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

		/* Update number of iterations */
		limits.iters_current += 1;
	}

	return best;
}
