/*
 * NEEMP - eem.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013, 2014
 *
 * */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

#ifdef USE_MKL
#include <mkl.h>
#endif /* USE_MKL */

#include "eem.h"
#include "neemp.h"
#include "settings.h"
#include "subset.h"
#include "structures.h"

extern const struct training_set ts;
extern const struct settings s;

#ifdef USE_MKL
static void fill_EEM_matrix_packed(float * const A, const struct molecule * const m, const struct kappa_data * const kd);
#else
static void fill_EEM_matrix_full(float * const A, const struct molecule * const m, const struct kappa_data * const kd);
static int GEM_solver(float * const A, float * const b, int n);
#endif /* USE_MKL */

#ifdef NOT_USED
static void print_matrix_packed(const float * const A, int n);
static void print_matrix_full(const float * const A, int n);

/* Print matrix in packed storage format */
static void print_matrix_packed(const float * const A, int n) {

	#define U_IDX(x, y) (x + (y * (y + 1))/2)
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < i; j++)
			printf("       ");
		for(int j = i; j < n; j++)
			printf("%6.4f ", A[U_IDX(i, j)]);

		printf("\n");
	}
	#undef U_IDX
}

/* Print matrix in packed storage format */
static void print_matrix_full(const float * const A, int n) {

	#define IDX(x, y) (x * n + y)
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++)
			fprintf(stderr, "%6.3f ", A[IDX(i, j)]);

		fprintf(stderr, "\n");
	}
	#undef IDX
}

#endif /* NOT_USED */

#ifdef USE_MKL
/* Use packed storage scheme to fill EEM matrix */
static void fill_EEM_matrix_packed(float * const A, const struct molecule * const m, const struct kappa_data * const kd) {

	assert(A != NULL);
	assert(m != NULL);
	assert(kd != NULL);

	const int n = m->atoms_count;

	/* Following #define works only for i <= j */
	#define U_IDX(x, y) (x + (y * (y + 1))/2)

	/* Fill the upper half of the n * n block */
	for(int i = 0; i < n; i++) {
		A[U_IDX(i, i)] = kd->parameters_beta[get_atom_type_idx(&m->atoms[i])];
		for(int j = i + 1; j < n; j++) {
			if(s.mode == MODE_PARAMS)
				A[U_IDX(i, j)] = (float) (kd->kappa * m->atoms[i].rdists[j]);
			else
				A[U_IDX(i, j)] = (float) (kd->kappa * rdist(&m->atoms[i], &m->atoms[j]));
		}
	}

	/* Fill last column */
	for(int i = 0; i < n; i++)
		A[U_IDX(i, n)] = 1.0f;

	/* Set the bottom right element to zero */
	A[U_IDX(n, n)] = 0.0f;
	#undef U_IDX
}
#else
/* Fill the whole standard EEM matrix */
static void fill_EEM_matrix_full(float * const A, const struct molecule * const m, const struct kappa_data * const kd) {

	assert(A != NULL);
	assert(m != NULL);
	assert(kd != NULL);

	const int n = m->atoms_count;

	#define IDX(x, y) (x * (n + 1) + y)
	/* Fill the n x n block */
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++)
			if(i == j)
				A[IDX(i, i)] = kd->parameters_beta[get_atom_type_idx(&m->atoms[i])];
			else {
				if(s.mode == MODE_PARAMS)
					A[IDX(i, j)] = (float) (kd->kappa * m->atoms[i].rdists[j]);
				else
					A[IDX(i, j)] = (float) (kd->kappa * rdist(&m->atoms[i], &m->atoms[j]));
			}
	}

	/* Fill last column */
	for(int i = 0; i < n; i++)
		A[IDX(i, n)] = 1.0f;

	/* Fill last row */
	for(int i = 0; i < n; i++)
		A[IDX(n, i)] = 1.0f;

	/* Set the bottom right element to zero */
	A[IDX(n, n)] = 0.0f;
	#undef IDX
}

/* Solve the system of linear eqns using GEM w/ partial pivoting */
static int GEM_solver(float * const A, float * const b, int n) {

	assert(A != NULL);
	assert(b != NULL);

	#define IDX(x, y) (x * n + y)
	/* GEM with partial pivoting */
	for(int i = 0; i < n; i++) {
		/* Find the largest pivot element */
		int best_row = i;
		for(int j = i + 1; j < n; j++)
			if(fabsf(A[IDX(best_row, i)]) < fabsf(A[IDX(j, i)]))
				best_row = j;

		/* Swap rows i and best_row */
		for(int j = 0; j < n; j++) {
			float tmp = A[IDX(best_row, j)];
			A[IDX(best_row, j)] = A[IDX(i, j)];
			A[IDX(i, j)] = tmp;
		}

		/* Swap elements in vector b in the same way */
		float tmp = b[best_row];
		b[best_row] = b[i];
		b[i] = tmp;

		/* Is matrix singular? */
		if(fabsf(A[IDX(i, i)]) < EPS)
			return i + 1;

		/* Actual elimination */
		for(int j = i + 1; j < n; j++) {
			const float mult = A[IDX(j, i)] / A[IDX(i, i)];
			for(int k = i + 1; k < n; k++)
				A[IDX(j, k)] -= mult * A[IDX(i, k)];

			b[j] -= mult * b[i];
		}
	}

	/* Backward substitution */
	for(int i = n - 1; i >= 0; i--) {
		b[i] /= A[IDX(i, i)];
		for(int j = i - 1; j >= 0; j--)
			b[j] -= A[IDX(j, i)] * b[i];
	}
	#undef IDX

	return 0;
}
#endif /* USE_MKL */

/* Calculate charges for a particular kappa_data structure */
void calculate_charges(struct subset * const ss, struct kappa_data * const kd) {

	assert(ss != NULL);
	assert(kd != NULL);

	/* Compute starting indices for storing the charges of each molecule. These are
	 * needed to guarantee the independence of the for loop iterations */
	int starts[ts.molecules_count];
	starts[0] = 0;
	for(int i = 1; i < ts.molecules_count; i++)
		starts[i] = starts[i - 1] + ts.molecules[i - 1].atoms_count;

	#pragma omp parallel for num_threads(ts.molecules_count < s.max_threads ? ts.molecules_count : s.max_threads)
	for(int i = 0; i < ts.molecules_count; i++) {
		#define MOLECULE ts.molecules[i]
		const int n = MOLECULE.atoms_count;

		#ifdef USE_MKL
		float *A = (float *) mkl_malloc(((n + 1) * (n + 2)) / 2 * sizeof(float), 32);
		float *b = (float *) mkl_malloc((n + 1) * sizeof(float), 32);
		#else
		float *A = (float *) malloc((n + 1) * (n + 1) * sizeof(float));
		float *b = (float *) malloc((n + 1) * sizeof(float));
		#endif /* USE_MKL */

		if(!A || !b)
			EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for EEM system.\n");

		#ifdef USE_MKL
		fill_EEM_matrix_packed(A, &ts.molecules[i], kd);
		#else
		fill_EEM_matrix_full(A, &ts.molecules[i], kd);
		#endif /* USE_MKL */

		/* Fill vector b */
		for(int j = 0; j < n; j++)
			b[j] = - kd->parameters_alpha[get_atom_type_idx(&MOLECULE.atoms[j])];

		b[n] = MOLECULE.sum_of_charges;

		/* Solve EEM system */
		#ifdef USE_MKL

		int ipiv[n + 1];

		/* int LAPACKE_sspsvx(int matrix_layout, char fact, char uplo, int n, int nrhs, const float *ap, float *afp,
		 * 	int *ipiv, const float *b, int ldb, float *x, int ldx, float *rcond, float *ferr, float *berr ); */

		float *Ap = (float *) mkl_malloc(((n + 1) * (n + 2)) / 2 * sizeof(float), 32);
		float *x = (float *) mkl_malloc((n + 1) * sizeof(float), 32);
		float rcond;
		float berr, ferr;

		LAPACKE_sspsvx(LAPACK_COL_MAJOR, 'N', 'U', n + 1, 1, A, Ap, ipiv, b, n + 1, x, n + 1, &rcond, &ferr, &berr);

		if(s.mode == MODE_CHARGES && rcond < WARN_MIN_RCOND)
			fprintf(stderr, "Ill-conditioned EEM system for molecule %s. Charges might be inaccurate.\n", MOLECULE.name);

		/* Store computed charges */
		for(int j = 0; j < n; j++)
			kd->charges[starts[i] + j] = x[j];

		#else

		GEM_solver(A, b, n + 1);

		/* Store computed charges */
		for(int j = 0; j < n; j++)
			kd->charges[starts[i] + j] = b[j];

		#endif /* USE_MKL */

		#undef MOLECULE

		/* Clean things up */
		#ifdef USE_MKL
		mkl_free(A);
		mkl_free(b);
		mkl_free(Ap);
		mkl_free(x);
		#else
		free(A);
		free(b);
		#endif /* USE_MKL */
	}
}
