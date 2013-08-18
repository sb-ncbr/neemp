/*
 * NEEMP - statistics.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "statistics.h"
#include "structures.h"
#include "subset.h"

extern const struct training_set ts;

/* Calculate R, RMSD and D for a kappa_data struct */
void calculate_statistics(struct subset * const ss, struct kappa_data * const kd) {

	assert(ss != NULL);
	assert(kd != NULL);

	int atoms_processed = 0;

	double D_sum_atoms = 0.0;
	double RMSD_sum_molecules = 0.0;
	double R_sum_molecules = 0.0;

	for(int i = 0; i < ts.molecules_count; i++) {
		#define MOLECULE ts.molecules[i]

		/* Calculate average charge from computed charges */
		double new_average_charge = 0.0;
		for(int j = 0; j < MOLECULE.atoms_count; j++)
			new_average_charge += kd->charges[atoms_processed + j];

		new_average_charge /= MOLECULE.atoms_count;

		double diff2 = 0.0;
		double sum1 = 0.0;
		double sum2 = 0.0;
		double sum3 = 0.0;

		for(int j = 0; j < MOLECULE.atoms_count; j++) {
			double x_xavg = MOLECULE.atoms[j].reference_charge - MOLECULE.average_charge;
			double y_yavg = kd->charges[atoms_processed + j] - new_average_charge;

			sum1 += x_xavg * y_yavg;
			sum2 += x_xavg * x_xavg;
			sum3 += y_yavg * y_yavg;

			double diff = kd->charges[atoms_processed + j] - MOLECULE.atoms[j].reference_charge;
			diff2 += diff * diff;
			D_sum_atoms += fabs(diff);
		}

		RMSD_sum_molecules += sqrt(diff2 / MOLECULE.atoms_count);
		R_sum_molecules += (sum1 * sum1) / (sum2 * sum3);

		atoms_processed += MOLECULE.atoms_count;
		#undef MOLECULE
	}

	kd->D = (float) (D_sum_atoms / ts.atoms_count);
	kd->R = (float) (R_sum_molecules / ts.molecules_count);
	kd->RMSD = (float) (RMSD_sum_molecules / ts.molecules_count);
}
