/*
 * NEEMP - subset.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013, 2014
 *
 * */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "bitarray.h"
#include "neemp.h"
#include "settings.h"
#include "structures.h"
#include "subset.h"

extern const struct settings s;
extern const struct training_set ts;

/* Initialize new subset structure from parent */
void ss_init(struct subset * const ss, const struct subset * const parent) {

	assert(ss != NULL);

	b_init(&ss->molecules, ts.molecules_count);
	ss->parent = parent;
	if(parent) {
		b_set_as(&ss->molecules, &parent->molecules);
	}
	else
		b_set_all(&ss->molecules);
}

/* Initialize empty subset with kappa_data */
void fill_ss(struct subset * const ss, int size) {

	assert(ss != NULL);

	ss->kappa_data_count = size;
	ss->data = (struct kappa_data*) calloc(ss->kappa_data_count, sizeof(struct kappa_data));
	if(!ss->data)
		EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for kappa data array.\n");
	for (int i = 0; i< ss->kappa_data_count; i++) {
		kd_init(&ss->data[i]);
	}
}

/* Allocate memory for the contents of kappa_data structure */
void kd_init(struct kappa_data * const kd) {

	assert(kd != NULL);

	kd->parameters_alpha = (float *) calloc(ts.atom_types_count, sizeof(float));
	kd->parameters_beta = (float *) calloc(ts.atom_types_count, sizeof(float));
	if(!kd->parameters_alpha || !kd->parameters_beta)
		EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for parameters array.\n");

	kd->charges = (float *) malloc(ts.atoms_count * sizeof(float));
	if(!kd->charges)
		EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for charges array.\n");

	kd->per_at_stats = (struct stats *) calloc(ts.atom_types_count, sizeof(struct stats));
	kd->per_molecule_stats = (struct stats *) calloc(ts.molecules_count, sizeof(struct stats));
}

/* Copy data from one kappa_data to another */
void kd_copy_parameters(struct kappa_data* from, struct kappa_data* to) {

	assert(from != NULL);
	assert(to != NULL);

	to->kappa = from->kappa;
	for (int i = 0; i < ts.atom_types_count; i++) {
		to->parameters_alpha[i] = from->parameters_alpha[i];
		to->parameters_beta[i] = from->parameters_beta[i];
	}
}

/* Copy statistics from one kappa_data to another */
void kd_copy_statistics(struct kappa_data* from, struct kappa_data* to) {

	assert(from != NULL);
	assert(to != NULL);

	to->full_stats.R = from->full_stats.R;
	to->full_stats.R2 = from->full_stats.R2;
	to->full_stats.R_w = from->full_stats.R_w;
	to->full_stats.RMSD = from->full_stats.RMSD;
	to->full_stats.spearman = from->full_stats.spearman;
	to->full_stats.D_avg = from->full_stats.D_avg;
	to->full_stats.D_max = from->full_stats.D_max;

	for (int i = 0; i < ts.atom_types_count; i++) {
		to->per_at_stats[i].R = from->per_at_stats[i].R;
		to->per_at_stats[i].R2 = from->per_at_stats[i].R2;
		to->per_at_stats[i].R_w = from->per_at_stats[i].R_w;
		to->per_at_stats[i].RMSD = from->per_at_stats[i].RMSD;
		to->per_at_stats[i].spearman = from->per_at_stats[i].spearman;
		to->per_at_stats[i].D_avg = from->per_at_stats[i].D_avg;
		to->per_at_stats[i].D_max = from->per_at_stats[i].D_max;
	}

}

/* Destroy contents of the kappa_data structure */
void kd_destroy(struct kappa_data * const kd) {

	assert(kd != NULL);

	free(kd->parameters_alpha);
	free(kd->parameters_beta);
	free(kd->charges);
	free(kd->per_at_stats);
	free(kd->per_molecule_stats);
}

/* Destroy contents of the subset */
void ss_destroy(struct subset * const ss) {

	assert(ss != NULL);

	for(int i = 0; i < ss->kappa_data_count; i++)
		kd_destroy(&ss->data[i]);

	b_destroy(&ss->molecules);
	free(ss->data);
}

/* Print loaded parameters */
void print_parameters(const struct kappa_data * const kd) {

	assert(kd != NULL);

	printf("\nLoaded parameters:\n");
	printf("K: %6.4f\n", kd->kappa);
	printf("Atom type             A       B\n");
	for(int i = 0; i < ts.atom_types_count; i++) {
		char buff[10];
		at_format_text(&ts.atom_types[i], buff);
		printf(" %s   \t%7.4f\t%7.4f\n", buff, kd->parameters_alpha[i], kd->parameters_beta[i]);
	}

	printf("\n");
}

/* Prints the parameters and associated stats */
void print_results(const struct subset * const ss) {

	assert(ss != NULL);
	assert(ss->best != NULL);

	printf("\nUsed molecules: %5d\n", b_count_bits(&ss->molecules));
	kd_print_stats(ss->best);

	printf("Atom type            A       B           R      R2      Sp      RMSD    D_avg   D_max\n");
	for(int i = 0; i < ts.atom_types_count; i++) {
		char buff[10];
		at_format_text(&ts.atom_types[i], buff);
		printf(" %s  \t%6.4f\t%6.4f      %6.4f  %6.4f  %6.4f    %6.4f   %6.4f  %6.4f\n",
			buff, ss->best->parameters_alpha[i], ss->best->parameters_beta[i],
			ss->best->per_at_stats[i].R, ss->best->per_at_stats[i].R2,
			ss->best->per_at_stats[i].spearman, ss->best->per_at_stats[i].RMSD,
			ss->best->per_at_stats[i].D_avg, ss->best->per_at_stats[i].D_max);
	}

	printf("\n");
}

/* Prints the parameters and associated stats */
void kd_print_results(const struct kappa_data * const kd) {

	assert(kd != NULL);

	kd_print_stats(kd);

	printf("Atom type            A       B           R      R2      Sp      RMSD    D_avg   D_max\n");
	for(int i = 0; i < ts.atom_types_count; i++) {
		char buff[10];
		at_format_text(&ts.atom_types[i], buff);
		printf(" %s  \t%6.4f\t%6.4f      %6.4f  %6.4f  %6.4f    %6.4f   %6.4f  %6.4f\n",
			buff, kd->parameters_alpha[i], kd->parameters_beta[i],
			kd->per_at_stats[i].R, kd->per_at_stats[i].R2,
			kd->per_at_stats[i].spearman, kd->per_at_stats[i].RMSD,
			kd->per_at_stats[i].D_avg, kd->per_at_stats[i].D_max);
	}

	printf("\n");
}

/* Return the value which we selected for sorting */
float kd_sort_by_return_value(const struct kappa_data * const kd) {

	assert(kd != NULL);

	switch(s.sort_by) {
		case SORT_R:
			return kd->full_stats.R;
		case SORT_R2:
			return kd->full_stats.R2;
		case SORT_RW:
			return kd->full_stats.R_w;
		case SORT_SPEARMAN:
			return kd->full_stats.spearman;
		case SORT_RMSD:
			return kd->full_stats.RMSD;
		case SORT_RMSD_AVG:
			return kd->full_stats.RMSD_avg;
		case SORT_D_AVG:
			return kd->full_stats.D_avg;
		case SORT_D_MAX:
			return kd->full_stats.D_max;
		default:
			/* Something bad happened */
			assert(0);
	}
}

/* Return the value which we selected for sorting from per atom type stats*/
float kd_sort_by_return_value_per_atom(const struct kappa_data * const kd, int i) {

	assert(kd != NULL);

	switch (s.sort_by) {
		case SORT_R:
			return kd->per_at_stats[i].R;
		case SORT_R2:
			return kd->per_at_stats[i].R2;
		case SORT_RW:
			return kd->per_at_stats[i].R_w;
		case SORT_SPEARMAN:
			return kd->per_at_stats[i].spearman;
		case SORT_RMSD:
			return kd->per_at_stats[i].RMSD;
		case SORT_D_AVG:
			return kd->per_at_stats[i].D_avg;
		case SORT_D_MAX:
			return kd->per_at_stats[i].D_max;
		default:
			assert(0);
	}
}

/* Determine if kd1 is better than kd2 in terms of the sort-by value */
int kd_sort_by_is_better(const struct kappa_data * const kd1, const struct kappa_data * const kd2) {

	assert(kd1 != NULL);
	assert(kd2 != NULL);

	/* Higher correlation is better, otherwise prefer lower value */
	if(s.sort_by == SORT_R || s.sort_by == SORT_R2 || s.sort_by == SORT_RW || s.sort_by == SORT_SPEARMAN)
		return kd_sort_by_return_value(kd1) > kd_sort_by_return_value(kd2);
	else
		return kd_sort_by_return_value(kd1) < kd_sort_by_return_value(kd2);
}

/* Determine if kd1 is much better or much worse in some element than kd2 in terms of the sort-by value per atom */
void kd_sort_by_is_much_better_per_atom(int* results_per_atom, const struct kappa_data * const kd1, const struct kappa_data * const kd2, float threshold) {

	assert(kd1 != NULL);
	assert(kd2 != NULL);

	for (int i = 0; i < ts.atom_types_count; i++) {
		if (kd_sort_by_return_value_per_atom(kd1, i)-kd_sort_by_return_value_per_atom(kd2, i) >= threshold)
			results_per_atom[i] = 1;
		else if (kd_sort_by_return_value_per_atom(kd2, i) - kd_sort_by_return_value_per_atom(kd1, i) >= threshold)
			results_per_atom[i] = -1;
		else
			results_per_atom[i] = 0;
	}
}

/* Determine if kd1 is better in some element than kd2 in terms of the sort-by value per atom */
int kd_sort_by_is_better_per_atom(const struct kappa_data * const kd1, const struct kappa_data * const kd2, int idx) {

	assert(kd1 != NULL);
	assert(kd2 != NULL);

	/* Higher correlation is better, otherwise prefer lower value */
	if(s.sort_by == SORT_R || s.sort_by == SORT_R2 || s.sort_by == SORT_RW || s.sort_by == SORT_SPEARMAN)
		return kd_sort_by_return_value_per_atom(kd1, idx) > kd_sort_by_return_value_per_atom(kd2, idx);
	else
		return kd_sort_by_return_value_per_atom(kd1, idx) < kd_sort_by_return_value_per_atom(kd2, idx);
}


/* Print all the statistics for the particular kappa data */
void kd_print_stats(const struct kappa_data * const kd) {

	assert(kd != NULL);

	/* Allocate buffer large enough to hold the entire message */
	char message[200];
	memset(message, 0, 200 * sizeof(char));

	snprintf(message, 200, "K: %6.4f |  R: %6.4f  R2: %6.4f  RW: %6.4f  Sp: %6.4f  RMSD: %6.4f  D_avg: %6.4f  D_max: %6.4f\n",
		kd->kappa, kd->full_stats.R, kd->full_stats.R2, kd->full_stats.R_w, kd->full_stats.spearman, kd->full_stats.RMSD, kd->full_stats.D_avg, kd->full_stats.D_max);

	printf("%s", message);
}
