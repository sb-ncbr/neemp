/*
 * NEEMP - statistics.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013, 2014
 *
 * */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "statistics.h"
#include "structures.h"
#include "subset.h"

extern const struct training_set ts;

/* Calculate R, RMSD and D, MSE and max/avg D per atom type for a given kappa_data struct */
void calculate_statistics(struct subset * const ss, struct kappa_data * const kd) {

	assert(ss != NULL);
	assert(kd != NULL);

	int atoms_processed = 0;

	double D_sum_atoms = 0.0;
	double MSE_sum_molecules = 0.0;
	double RMSD_sum_molecules = 0.0;
	double R_sum_molecules = 0.0;

	double *D_sum_atom_type = (double *) calloc(ts.atom_types_count, sizeof(double));

	/* Reset values left by previous iteration of the Brent's method */
	for(int i = 0; i < ts.atom_types_count; i++)
		kd->per_at_stats[i].D_max = 0.0f;

	for(int i = 0; i < ts.molecules_count; i++) {
		#define MOLECULE ts.molecules[i]

		/* Calculate average charge from computed charges */
		double new_average_charge = 0.0;
		for(int j = 0; j < MOLECULE.atoms_count; j++)
			new_average_charge += kd->charges[atoms_processed + j];

		new_average_charge /= MOLECULE.atoms_count;

		double diff2_sum = 0.0;
		double sum_xy = 0.0;
		double sum_xx = 0.0;
		double sum_yy = 0.0;

		double max_diff_per_molecule = 0.0;
		double D_sum_molecule = 0.0;

		for(int j = 0; j < MOLECULE.atoms_count; j++) {
			double x_xavg = MOLECULE.atoms[j].reference_charge - MOLECULE.average_charge;
			double y_yavg = kd->charges[atoms_processed + j] - new_average_charge;

			sum_xy += x_xavg * y_yavg;
			sum_xx += x_xavg * x_xavg;
			sum_yy += y_yavg * y_yavg;

			double diff = kd->charges[atoms_processed + j] - MOLECULE.atoms[j].reference_charge;
			diff2_sum += diff * diff;
			D_sum_atoms += fabs(diff);
			D_sum_molecule += fabs(diff);

			if(max_diff_per_molecule < fabs(diff))
				max_diff_per_molecule = fabs(diff);

			const int atom_type_idx = get_atom_type_idx(&MOLECULE.atoms[j]);

			if(fabsf(kd->per_at_stats[atom_type_idx].D_max) < fabs(diff))
				kd->per_at_stats[atom_type_idx].D_max = (float) fabs(diff);

			D_sum_atom_type[atom_type_idx] += fabs(diff);
		}

		R_sum_molecules += (sum_xy * sum_xy) / (sum_xx * sum_yy);
		RMSD_sum_molecules += sqrt(diff2_sum / MOLECULE.atoms_count);
		MSE_sum_molecules += diff2_sum;

		kd->per_molecule_stats[i].R = (float) ((sum_xy * sum_xy) / (sum_xx * sum_yy));
		kd->per_molecule_stats[i].RMSD = (float) sqrt(diff2_sum / MOLECULE.atoms_count);
		kd->per_molecule_stats[i].MSE = (float) diff2_sum / MOLECULE.atoms_count;
		kd->per_molecule_stats[i].D_avg = (float) (D_sum_molecule / MOLECULE.atoms_count);
		kd->per_molecule_stats[i].D_max = (float) max_diff_per_molecule;

		atoms_processed += MOLECULE.atoms_count;
		#undef MOLECULE
	}

	kd->full_stats.R = (float) (R_sum_molecules / ts.molecules_count);
	kd->full_stats.RMSD = (float) (RMSD_sum_molecules / ts.molecules_count);
	kd->full_stats.MSE = (float) (MSE_sum_molecules / ts.molecules_count);
	kd->full_stats.D_avg = (float) (D_sum_atoms / ts.atoms_count);

	/* Calculate per atom type statistics */
	for(int i = 0; i < ts.atom_types_count; i++)
		kd->per_at_stats[i].D_avg = (float) D_sum_atom_type[i] / ts.atom_types[i].atoms_count;

	kd->full_stats.D_max = kd->per_at_stats[0].D_max;
	for(int i = 1; i < ts.atom_types_count; i++)
		if(kd->full_stats.D_max < kd->per_at_stats[i].D_max)
			kd->full_stats.D_max = kd->per_at_stats[i].D_max;

	free(D_sum_atom_type);

	/* Compute starting indices for storing the charges of each molecule. These are
	 * needed to access individual charges. */
	int starts[ts.molecules_count];
	starts[0] = 0;
	for(int i = 1; i < ts.molecules_count; i++)
		starts[i] = starts[i - 1] + ts.molecules[i - 1].atoms_count;

	for(int i = 0; i < ts.atom_types_count; i++) {
		#define AT ts.atom_types[i]

		double avg_computed_charge_per_atom_type = 0.0;
		double avg_reference_charge_per_atom_type = 0.0;

		for(int j = 0; j < ts.atom_types[i].atoms_count; j++) {
			const int molecule_idx = AT.atoms_molecule_idx[j];
			const int atom_idx = AT.atoms_atom_idx[j];

			avg_computed_charge_per_atom_type += kd->charges[starts[molecule_idx] + atom_idx];
			avg_reference_charge_per_atom_type += ts.molecules[molecule_idx].atoms[atom_idx].reference_charge;
		}

		avg_computed_charge_per_atom_type /= AT.atoms_count;
		avg_reference_charge_per_atom_type /= AT.atoms_count;

		double diff2_sum = 0.0;
		double sum_xy = 0.0;
		double sum_xx = 0.0;
		double sum_yy = 0.0;

		for(int j  = 0; j < ts.atom_types[i].atoms_count; j++) {
			const int molecule_idx = AT.atoms_molecule_idx[j];
			const int atom_idx = AT.atoms_atom_idx[j];

			double diff = ts.molecules[molecule_idx].atoms[atom_idx].reference_charge - kd->charges[starts[molecule_idx] + atom_idx];
			diff2_sum += diff * diff;

			double x_xavg = ts.molecules[molecule_idx].atoms[atom_idx].reference_charge - avg_reference_charge_per_atom_type;
			double y_yavg = kd->charges[starts[molecule_idx] + atom_idx] - avg_computed_charge_per_atom_type;

			sum_xy += x_xavg * y_yavg;
			sum_xx += x_xavg * x_xavg;
			sum_yy += y_yavg * y_yavg;
		}

		kd->per_at_stats[i].R = (float) ((sum_xy * sum_xy) / (sum_xx * sum_yy));
		kd->per_at_stats[i].RMSD = (float) sqrt(diff2_sum / ts.atom_types[i].atoms_count);
		kd->per_at_stats[i].MSE = (float) (diff2_sum / ts.atom_types[i].atoms_count);
		#undef AT
	}
}

/* Check for abnormal charge differences */
void check_charges(const struct kappa_data * const kd) {

	assert(kd != NULL);

	int bad_molecules = 0;

	for(int i = 0; i < ts.molecules_count; i++) {
		if(kd->per_molecule_stats[i].R < WARN_MIN_R ||
		   kd->per_molecule_stats[i].RMSD > WARN_MAX_RMSD ||
		   kd->per_molecule_stats[i].MSE > WARN_MAX_MSE ||
		   kd->per_molecule_stats[i].D_avg > WARN_MAX_D_AVG ||
		   kd->per_molecule_stats[i].D_max > WARN_MAX_D_MAX) {
			fprintf(stderr, "Warning: Abnormal values for molecule %s\n", ts.molecules[i].name);
			fprintf(stderr, "R: %6.4f  RMSD: %4.2e  MSE: %4.2e  D_avg: %4.2e  D_max: %4.2e\n\n",
				kd->per_molecule_stats[i].R, kd->per_molecule_stats[i].RMSD, kd->per_molecule_stats[i].MSE,
				kd->per_molecule_stats[i].D_avg, kd->per_molecule_stats[i].D_max);

			bad_molecules++;
		}

	}
	fprintf(stderr, "Check charges: %d molecules with abnormal statistics found.\n", bad_molecules);
}