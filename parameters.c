/*
 * NEEMP - parameters.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#include <assert.h>
#include <mkl.h>

#include <stdlib.h>

#include "parameters.h"
#include "structures.h"

extern struct training_set ts;

/* Calculate parameters for give kappa_data structure */
void calculate_parameters(struct subset * const ss, struct kappa_data * const kd) {

	assert(ss != NULL);
	assert(kd != NULL);

	for(int i = 0; i < ts.atom_types_count; i++) {
		#define AT ts.atom_types[i]
		double *A = (double *) mkl_malloc(2 * AT.atoms_count * sizeof(double), 32);
		double *b = (double *) mkl_malloc(AT.atoms_count * sizeof(double), 32);

		for(int j = 0; j < AT.atoms_count; j++) {
			#define MOLECULE ts.molecules[AT.atoms_molecule_idx[j]]
			#define ATOM ts.molecules[AT.atoms_molecule_idx[j]].atoms[AT.atoms_atom_idx[j]]
			A[2 * j] = 1.0;
			A[2 * j + 1] = ATOM.reference_charge;
			b[j] = MOLECULE.electronegativity - kd->kappa * ATOM.y;
			#undef ATOM
			#undef MOLECULE
		}

		LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', AT.atoms_count, 2, 1, A, 2, b, 1);

		kd->parameters_alpha[i] = (float) b[0];
		kd->parameters_beta[i] = (float) b[1];
		#undef AT

		mkl_free(A);
		mkl_free(b);
	}
}
