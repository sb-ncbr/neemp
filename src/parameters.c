/*
 * NEEMP - parameters.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013, 2014
 *
 * */

#include <assert.h>
#include <stdlib.h>

#ifdef USE_MKL
#include <mkl.h>
#endif /* USE_MKL */

#include "neemp.h"
#include "parameters.h"
#include "structures.h"

extern struct training_set ts;

static inline int is_molecule_enabled(const struct subset * const ss, int mol_idx);

static inline int is_molecule_enabled(const struct subset * const ss, int mol_idx) {

	assert(ss != NULL);

	return b_get(&ss->molecules, mol_idx);
}

#ifdef USE_MKL
/* Calculate parameters for give kappa_data structure */
void calculate_parameters(struct subset * const ss, struct kappa_data * const kd) {

	assert(ss != NULL);
	assert(kd != NULL);

	for(int i = 0; i < ts.atom_types_count; i++) {
		#define AT ts.atom_types[i]

		double *A = (double *) mkl_malloc(2 * AT.atoms_count * sizeof(double), 32);
		/* b has to have at least 2 elements since b[1] is used to store beta parameter.
		 * If atoms_count is 1 (which is unlikely but possible) we would probably crash... */
		double *b = (double *) mkl_malloc((1 + AT.atoms_count) * sizeof(double), 32);

		/* If some molecule is disabled, skip its atoms.
		 * To avoid having empty (zero) rows, change the index so that all used atoms are
		 * at the beginning. */
		int disabled_count = 0;
		for(int j = 0; j < AT.atoms_count;j++) {
			#define MOLECULE ts.molecules[AT.atoms_molecule_idx[j]]
			#define ATOM ts.molecules[AT.atoms_molecule_idx[j]].atoms[AT.atoms_atom_idx[j]]

			if(is_molecule_enabled(ss, AT.atoms_molecule_idx[j])) {
				A[2 * (j - disabled_count)] = 1.0;
				A[2 * (j - disabled_count) + 1] = ATOM.reference_charge;
				b[j - disabled_count] = MOLECULE.electronegativity - kd->kappa * ATOM.y;
			} else
				disabled_count++;

			#undef ATOM
			#undef MOLECULE
		}

		int info = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', AT.atoms_count - disabled_count, 2, 1, A, 2, b, 1);
		if(info)
			EXIT_ERROR(RUN_ERROR, "%s", "The least squares method failed.\n");

		kd->parameters_alpha[i] = (float) b[0];
		kd->parameters_beta[i] = (float) b[1];
		#undef AT

		mkl_free(A);
		mkl_free(b);
	}
}
#else
void calculate_parameters(struct subset * const ss, struct kappa_data * const kd) {

	assert(ss != NULL);
	assert(kd != NULL);

	for(int i = 0; i < ts.atom_types_count; i++) {
		#define AT ts.atom_types[i]

		double sum1 = 0.0;
		double sum2 = 0.0;
		double sum3 = 0.0;
		double sum4 = 0.0;

		int disabled_count = 0;
		for(int j = 0; j < AT.atoms_count; j++) {
			#define MOLECULE ts.molecules[AT.atoms_molecule_idx[j]]
			#define ATOM ts.molecules[AT.atoms_molecule_idx[j]].atoms[AT.atoms_atom_idx[j]]
			if(!is_molecule_enabled(ss, AT.atoms_molecule_idx[j])) {
				disabled_count++;
				continue;
			}

			double x = MOLECULE.electronegativity - kd->kappa * ATOM.y;

			sum1 += ATOM.reference_charge * x;
			sum2 += ATOM.reference_charge;
			sum3 += x;
			sum4 += ATOM.reference_charge * ATOM.reference_charge;

			#undef ATOM
			#undef MOLECULE
		}

		const int atoms_used = AT.atoms_count - disabled_count;
		kd->parameters_beta[i] = (float) ((atoms_used * sum1  - sum2 * sum3) / (atoms_used * sum4 - sum2 * sum2));
		kd->parameters_alpha[i] = (float) ((sum3 - kd->parameters_beta[i] * sum2) / atoms_used);
		#undef AT
	}
}
#endif /* USE_MKL */
