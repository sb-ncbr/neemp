/*
 * NEEMP - kappa.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "eem.h"
#include "kappa.h"
#include "neemp.h"
#include "parameters.h"
#include "settings.h"
#include "subset.h"
#include "statistics.h"
#include "structures.h"

extern const struct training_set ts;
extern const struct settings s;

static void full_scan(struct subset * const ss);
static void brent(struct subset * const ss);
static void perform_calculations(struct subset * const ss, struct kappa_data * const kd);

/* Perform all three steps for one value of kappa */
static void perform_calculations(struct subset * const ss, struct kappa_data * const kd) {

	assert(ss != NULL);
	assert(kd != NULL);

	calculate_parameters(ss, kd);
	calculate_charges(ss, kd);
	calculate_statistics(ss, kd);
}

/* Perform full scan */
static void full_scan(struct subset * const ss) {

	assert(ss != NULL);

	for(int i = 0; i < ss->kappa_data_count - 1; i++) {
		ss->data[i].kappa = i * s.full_scan_precision;
		perform_calculations(ss, &ss->data[i]);

		if(s.verbosity >= VERBOSE_KAPPA) {
			printf("F> ");
			kd_print_stats(&ss->data[i]);
		}
	}
}

/* Run full scan followed by the Brent's method to polish the result;
 * note that it can be only used for R */
static void brent(struct subset * const ss) {

	assert(ss != NULL);

	full_scan(ss);

	if(ss->kappa_data_count < 3)
		EXIT_ERROR(RUN_ERROR, "%s", "Cannot determine the initial inverval for the Brent's method.\n");

	/* Find the best so far */
	int best_idx = 0;
	for(int i = 1; i < ss->kappa_data_count - 1; i++)
		if(ss->data[i].stats.R > ss->data[best_idx].stats.R)
			best_idx = i;

	/* Create initial inverval for Brent */
	int left_idx = -1;
	int right_idx = -1;

	/* The left point is the best */
	if(best_idx == 0) {
		left_idx = 0;
		right_idx = 1;
	} else if(best_idx == ss->kappa_data_count - 2) {
		left_idx = ss->kappa_data_count - 3;
		right_idx = ss->kappa_data_count - 2;
	}
	else {
		left_idx = best_idx - 1;
		right_idx = best_idx + 1;
	}

	#define CGOLD 0.3819660F
	#define ZEPS 1.0e-10F
	#define SHFT(a, b, c, d) (a) = (b); (b) = (c); (c) = (d);
	#define BRENT_MAX_ITERS 20
	#define SIGN(a, b) ((b) >= 0.0f ? fabsf(a) : -fabsf(a))
	#define TOLERANCE 0.001F

	#define KAPPA_DATA_BRENT ss->data[ss->kappa_data_count - 1]

	/* Initial interval values */
	float a = ss->data[left_idx].kappa;
	float b = ss->data[right_idx].kappa;

	/* Abscissa */
	float x, w, v, u;
	/* Functional values */
	float fx, fw, fv, fu;
	/* Trial parabolic fit */
	float p, q, r;
	/* Accepting tolerances */
	float tol1, tol2;

	float d = 0.0f;
	float e = 0.0f;

	x = w = v = 0.5f * (a + b);

	KAPPA_DATA_BRENT.kappa = x;
	perform_calculations(ss, &KAPPA_DATA_BRENT);
	/* The code is for minimization, so take the negative of R */
	fw = fv = fx=  - KAPPA_DATA_BRENT.stats.R;

	if(s.verbosity >= VERBOSE_KAPPA) {
			printf("B> ");
			kd_print_stats(&KAPPA_DATA_BRENT);
	}

	for(int i = 0; i < BRENT_MAX_ITERS; i++) {
		float xm = 0.5f * (a + b);
		tol2 = 2.0f * (tol1 = TOLERANCE * fabsf(x) + ZEPS);

		/* Check for acceptance of the solution */
		if(fabsf(x - xm) <= (tol2 - 0.5f * (b - a)))
			break;

		/* Try parabolic fit */
		if (fabsf(e) > tol1) {
			r = (x - w) * (fx - fv);
			q = (x - v) * (fx - fw);
			p = (x - v) * q - (x - w) * r;
			q = 2.0f * (q - r);
			if(q > 0.0f)
				p = -p;

			float etemp = e;
			q = fabsf(q);
			e = d;

			if (fabsf(p) >= fabsf(0.5f * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
				d = CGOLD * (e = (x >= xm ? a - x : b - x));
			else {
				d = p / q;
				u = x + d;
				if(u - a < tol2 || b - u < tol2)
					d = SIGN(tol1, xm - x);
			}
		} else
			d = CGOLD * (e = (x >= xm ? a - x : b - x));

		u = (fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
		KAPPA_DATA_BRENT.kappa = u;
		perform_calculations(ss, &KAPPA_DATA_BRENT);
		/* The code is for minimization, so take the negative of R */
		fu = - KAPPA_DATA_BRENT.stats.R;

		if(s.verbosity >= VERBOSE_KAPPA) {
				printf("B> ");
				kd_print_stats(&KAPPA_DATA_BRENT);
		}

		/* Housekeeping */
		if(fu <= fx) {
			if(u >= x)
				a = x;
			else
				b = x;

			SHFT(v, w, x, u)
			SHFT(fv, fw, fx, fu)
		} else {
			if(u < x)
				a = u;
			else
				b = u;

			if(fu <= fw || fabsf(w - x) < TOLERANCE) {
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			} else if(fu <= fv || fabsf(v - x) < TOLERANCE || fabsf(v - w) < TOLERANCE) {
				v = u;
				fv = fu;
			}
		}
	}
	#undef KAPPA_DATA_BRENT
	#undef BRENT_MAX_ITERS
	#undef CGOLD
	#undef ZEPS
	#undef SHFT
}

/* Finds the best parameters for a particular subset */
void find_the_best_parameters_for_subset(struct subset * const ss) {

	assert(ss != NULL);

	if(s.kappa_set > 1e-10) {
		ss->kappa_data_count = 1;
		ss->data = (struct kappa_data *) malloc(1 * sizeof(struct kappa_data));
		kd_init(&ss->data[0]);
		ss->data[0].kappa = s.kappa_set;

		perform_calculations(ss, &ss->data[0]);
		ss->best = &ss->data[0];
	}
	else {
		/* The last item in ss->data array is reserved for Brent whether it's used or not */
		ss->kappa_data_count = 1 + (int) (s.kappa_max / s.full_scan_precision);
		ss->data = (struct kappa_data *) calloc(ss->kappa_data_count, sizeof(struct kappa_data));
		if(!ss->data)
			EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for kappa data array.\n");

		for(int i = 0; i < ss->kappa_data_count; i++)
			kd_init(&ss->data[i]);

		if(s.full_scan_only)
			full_scan(ss);
		else
			brent(ss);

		/* Determine the best parameters for computed data */
		ss->best = &ss->data[0];

		if(s.full_scan_only) {
			for(int i = 0; i < ss->kappa_data_count - 1; i++)
				if(kd_sort_by_is_better(&ss->data[i], ss->best))
					ss->best = &ss->data[i];
		}
		else {
			/* If Brent is used, the maximum is stored in the last item */
			ss->best = &ss->data[ss->kappa_data_count - 1];
		}

	}
}
