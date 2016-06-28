/* Copyright 2013-2016 Tomas Racek (tom@krab1k.net)
 *
 * This file is part of NEEMP.
 *
 * NEEMP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * NEEMP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with NEEMP. If not, see <http://www.gnu.org/licenses/>.
 */

#include <assert.h>
#include <stdlib.h>

#ifdef USE_MKL
#include <mkl.h>
#else
extern void dgels_(const char* trans, const int *m, const int *n, const int *nrhs, double *a, const int *lda, double *b, const int *ldb, double *work, const int *lwork, int *info);
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

/* Calculate parameters for give kappa_data structure */
void calculate_parameters(struct subset * const ss, struct kappa_data * const kd) {

	assert(ss != NULL);
	assert(kd != NULL);

	for(int i = 0; i < ts.atom_types_count; i++) {
		#define AT ts.atom_types[i]

		#ifdef USE_MKL
		double *A = (double *) mkl_malloc(2 * AT.atoms_count * sizeof(double), 64);
		/* b has to have at least 2 elements since b[1] is used to store beta parameter.
		 * If atoms_count is 1 (which is unlikely but possible) we would probably crash... */
		double *b = (double *) mkl_malloc((1 + AT.atoms_count) * sizeof(double), 64);
		#else
		double *A = (double *) malloc(2 * AT.atoms_count * sizeof(double));
		double *b = (double *) malloc((1 + AT.atoms_count) * sizeof(double));
		#endif /* USE_MKL */

		/* If some molecule is disabled, skip its atoms.
		 * To avoid having empty (zero) rows, change the index so that all used atoms are
		 * at the beginning. */
		int disabled_count = 0;
		for(int j = 0; j < AT.atoms_count;j++) {
			#define MOLECULE ts.molecules[AT.atoms_molecule_idx[j]]
			#define ATOM ts.molecules[AT.atoms_molecule_idx[j]].atoms[AT.atoms_atom_idx[j]]

			if(is_molecule_enabled(ss, AT.atoms_molecule_idx[j])) {
				A[j - disabled_count] = 1.0;
				b[j - disabled_count] = MOLECULE.electronegativity - kd->kappa * ATOM.y;
			} else
				disabled_count++;

			#undef ATOM
			#undef MOLECULE
		}

		int border = AT.atoms_count - disabled_count;
		disabled_count = 0;

		for(int j = 0; j < AT.atoms_count;j++) {
			#define MOLECULE ts.molecules[AT.atoms_molecule_idx[j]]
			#define ATOM ts.molecules[AT.atoms_molecule_idx[j]].atoms[AT.atoms_atom_idx[j]]
			if(is_molecule_enabled(ss, AT.atoms_molecule_idx[j]))
				A[border + j - disabled_count] = ATOM.reference_charge;
			else
				disabled_count++;
			#undef ATOM
			#undef MOLECULE
		}

		char trans = 'N';
		int m = AT.atoms_count - disabled_count;
		int n = 2;
		int nrhs = 1;
		int lda = m;
		int ldb = m;
		int info;

		#ifdef USE_MKL
		info = LAPACKE_dgels(LAPACK_COL_MAJOR, trans, m, n, nrhs, A, lda, b, ldb);
		#else
		int lwork = -1;
		double work_opt;

		/* Check for optimal value of lwork */
		dgels_(&trans, &m, &n, &nrhs, A, &lda, b, &ldb, &work_opt, &lwork, &info);
		lwork = (int) work_opt;
		double work[lwork];
		dgels_(&trans, &m, &n, &nrhs, A, &lda, b, &ldb, work, &lwork, &info);
		#endif /* USE_MKL */

		if(info)
			EXIT_ERROR(RUN_ERROR, "%s", "The least squares method failed.\n");

		kd->parameters_alpha[i] = (float) b[0];
		kd->parameters_beta[i] = (float) b[1];
		#undef AT

		#ifdef USE_MKL
		mkl_free(A);
		mkl_free(b);
		#else
		free(A);
		free(b);
		#endif /* USE_MKL */
	}
}
