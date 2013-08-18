/*
 * NEEMP - eem.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#include <assert.h>
#include <mkl.h>
#include <stdlib.h>

#include "eem.h"
#include "neemp.h"
#include "subset.h"
#include "structures.h"

extern const struct training_set ts;

static void fill_EEM_matrix_packed(float * const A, const struct molecule * const m, const struct kappa_data * const kd);

/* Use packed storage scheme to fill EEM matrix */
static void fill_EEM_matrix_packed(float * const A, const struct molecule * const m, const struct kappa_data * const kd) {

	assert(A != NULL);
	assert(m != NULL);
	assert(kd != NULL);

	const int n = m->atoms_count;

	/* Following #define works only for i <= j */
	#define U_IDX(x, y) (x + (y * (y + 1))/2)

	/* Fill upper half of the n * n block */
	for(int i = 0; i < n; i++) {
		A[U_IDX(i, i)] = kd->parameters_beta[get_atom_type_idx(&m->atoms[i])];
		for(int j = i + 1; j < n; j++) {
			A[U_IDX(i, j)] = (float) (kd->kappa * m->atoms[i].rdists[j]);
		}
	}

	/* Fill last column */
	for(int i = 0; i < n; i++)
		A[U_IDX(i, n)] = 1.0f;

	/* Set the bottom right element to zero */
	A[U_IDX(n, n)] = 0.0f;

	#undef U_IDX
}

/* Calculate charges for a particular kappa_data structure */
void calculate_charges(struct subset * const ss, struct kappa_data * const kd) {

	assert(ss != NULL);
	assert(kd != NULL);

	int atoms_processed = 0;

	for(int i = 0; i < ts.molecules_count; i++) {
		#define MOLECULE ts.molecules[i]
		const int n = MOLECULE.atoms_count;
		float *A = (float *) mkl_malloc(((n + 1) * (n + 2)) / 2 * sizeof(float), 32);
		float *b = (float *) mkl_malloc((n + 1) * sizeof(float), 32);
		if(!A || !b)
			EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for EEM system.\n");

		fill_EEM_matrix_packed(A, &ts.molecules[i], kd);

		/* Fill vector b */
		for(int j = 0; j < n; j++)
			b[j] = -1 * kd->parameters_alpha[get_atom_type_idx(&MOLECULE.atoms[j])];

		b[n] = MOLECULE.sum_of_charges;

		int ipiv[n + 1];

		LAPACKE_sspsv(LAPACK_COL_MAJOR, 'U', n + 1, 1, A, ipiv, b, n + 1);

		/* Store computed charges */
		for(int j = 0; j < n; j++)
			kd->charges[atoms_processed + j] = b[j];

		atoms_processed += n;

		mkl_free(A);
		mkl_free(b);
		#undef MOLECULE
	}
}
