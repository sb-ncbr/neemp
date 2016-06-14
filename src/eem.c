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
#else
extern void sspsvx_(char *fact, char *uplo, int *n, int *nrhs, const float *ap, float *afp, int *ipiv, const float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr, float *berr, float *work, int *iwork, int *info);
#endif /* USE_MKL */

#include "eem.h"
#include "neemp.h"
#include "settings.h"
#include "subset.h"
#include "structures.h"

extern const struct training_set ts;
extern const struct settings s;


static void fill_EEM_matrix_packed(float * const A, const struct molecule * const m, const struct kappa_data * const kd);

#ifdef NOT_USED
static void print_matrix_packed(const float * const A, long int n);
static void print_matrix_full(const float * const A, long int n);

/* Print matrix in packed storage format */
static void print_matrix_packed(const float * const A, long int n) {

	assert(A != NULL);

	#define U_IDX(x, y) (x + (y * (y + 1))/2)
	for(long int i = 0; i < n; i++) {
		for(long int j = 0; j < i; j++)
			printf("       ");
		for(long int j = i; j < n; j++)
			printf("%6.4f ", A[U_IDX(i, j)]);

		printf("\n");
	}
	#undef U_IDX
}

/* Print matrix in packed storage format */
static void print_matrix_full(const float * const A, int n) {

	assert(A != NULL);

	#define IDX(x, y) (x * n + y)
	for(long int i = 0; i < n; i++) {
		for(long int j = 0; j < n; j++)
			fprintf(stderr, "%6.3f ", A[IDX(i, j)]);

		fprintf(stderr, "\n");
	}
	#undef IDX
}

#endif /* NOT_USED */

/* Use packed storage scheme to fill EEM matrix */
static void fill_EEM_matrix_packed(float * const A, const struct molecule * const m, const struct kappa_data * const kd) {

	assert(A != NULL);
	assert(m != NULL);
	assert(kd != NULL);

	const long int n = m->atoms_count;

	/* Following #define works only for i <= j */
	#define U_IDX(x, y) (x + (y * (y + 1))/2)

	/* Fill the upper half of the n * n block */
	for(long int i = 0; i < n; i++) {
		A[U_IDX(i, i)] = kd->parameters_beta[get_atom_type_idx(&m->atoms[i])];
		for(long int j = i + 1; j < n; j++) {
			if(s.mode == MODE_PARAMS)
				A[U_IDX(i, j)] = (float) (kd->kappa * m->atoms[i].rdists[j]);
			else
				A[U_IDX(i, j)] = (float) (kd->kappa * rdist(&m->atoms[i], &m->atoms[j]));
		}
	}

	/* Fill last column */
	for(long int i = 0; i < n; i++)
		A[U_IDX(i, n)] = 1.0f;

	/* Set the bottom right element to zero */
	A[U_IDX(n, n)] = 0.0f;
	#undef U_IDX
}

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
	int nt = s.max_threads;
	if (s.params_method == PARAMS_DE)
		nt /= s.de_threads;
	if (s.params_method == PARAMS_GM)
		nt /= s.gm_threads;

	int nthreads = ts.molecules_count < nt ? ts.molecules_count : nt;
	#pragma omp parallel for num_threads(nthreads)
	for(int i = 0; i < ts.molecules_count; i++) {
		#define MOLECULE ts.molecules[i]
		const int n = MOLECULE.atoms_count;

		#ifdef USE_MKL
		float *A = (float *) mkl_malloc(((n + 1) * (n + 2)) / 2 * sizeof(float), 64);
		float *b = (float *) mkl_malloc((n + 1) * sizeof(float), 64);
		float *Ap = (float *) mkl_malloc(((n + 1) * (n + 2)) / 2 * sizeof(float), 64);
		float *x = (float *) mkl_malloc((n + 1) * sizeof(float), 64);
		#else
		float *A = (float *) malloc((n + 1) * (n + 2) / 2 * sizeof(float));
		float *b = (float *) malloc((n + 1) * sizeof(float));
		float *Ap = (float *) malloc(((n + 1) * (n + 2)) / 2 * sizeof(float));
		float *x = (float *) malloc((n + 1) * sizeof(float));
		#endif /* USE_MKL */

		if(!A || !b || !Ap || !x)
			EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for EEM system.\n");

		fill_EEM_matrix_packed(A, &ts.molecules[i], kd);

		/* Fill vector b */
		for(int j = 0; j < n; j++)
			b[j] = - kd->parameters_alpha[get_atom_type_idx(&MOLECULE.atoms[j])];

		b[n] = MOLECULE.sum_of_charges;

		/* Solve EEM system */

		int ipiv[n + 1];
		float rcond;
		float berr, ferr;
		char fact = 'N';
		char uplo = 'U';
		int nn = n + 1;
		int nrhs = 1;
		int ldb = n + 1;
		int ldx = n + 1;

		#ifdef USE_MKL
		LAPACKE_sspsvx(LAPACK_COL_MAJOR, fact, uplo, nn, nrhs, A, Ap, ipiv, b, ldb, x, ldx, &rcond, &ferr, &berr);
		#else
		int info;
		int iwork[n + 1];
		float work[3 * (n + 1)];

		sspsvx_(&fact, &uplo, &nn, &nrhs, A, Ap, ipiv, b, &ldb, x, &ldx, &rcond, &ferr, &berr, work, iwork, &info);
		#endif /* USE_MKL */

		kd->per_molecule_stats[i].cond = 1 / rcond;

		if(s.mode == MODE_CHARGES && (1 / rcond) > WARN_MAX_COND)
			fprintf(stderr, "Ill-conditioned EEM system for molecule %s. Charges might be inaccurate.\n", MOLECULE.name);

		/* Store computed charges */
		for(int j = 0; j < n; j++)
			kd->charges[starts[i] + j] = x[j];

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
		free(Ap);
		free(x);
		#endif /* USE_MKL */
	}
}
