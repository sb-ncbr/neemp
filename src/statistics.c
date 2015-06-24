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
#include "neemp.h"
#include "settings.h"
#include "statistics.h"
#include "structures.h"
#include "subset.h"

extern const struct training_set ts;
extern const struct settings s;

/* Adjust weights for atom types if not specified by user */
void adjust_weights(struct subset * const ss) {

	for(int i = 0; i < ts.atom_types_count; i++) {
		ss->weights[i] = ((float) ts.atoms_count) / (ts.atom_types[i].atoms_count * ts.atom_types_count);
	}
}

static void set_total_R(struct kappa_data * const kd);
static void set_total_RMSD(struct kappa_data * const kd);
static void set_total_D_avg(struct kappa_data * const kd);
static void set_total_D_max(struct kappa_data * const kd);
static void set_total_R_weighted(struct subset * const ss, struct kappa_data * const kd);

static void set_per_at_R(struct kappa_data * const kd);
static void set_per_at_RMSD(struct kappa_data * const kd);
static void set_per_at_D_avg(struct kappa_data * const kd);
static void set_per_at_D_max(struct kappa_data * const kd);


/* Set total Pearson correlation coeff. squared for the whole set */
static void set_total_R(struct kappa_data * const kd) {

	assert(kd != NULL);

	int atoms_processed = 0;

	int bad_molecules = 0;
	double R_sum_molecules = 0.0;

	for(int i = 0; i < ts.molecules_count; i++) {
		#define MOLECULE ts.molecules[i]

		double average_calculated_charge = 0.0;
		for(int j = 0; j < MOLECULE.atoms_count; j++) {
			average_calculated_charge += kd->charges[atoms_processed + j];
		}

		average_calculated_charge /= MOLECULE.atoms_count;

		double cov_xy = 0.0;
		double cov_xx = 0.0;
		double cov_yy = 0.0;

		for(int j = 0; j < MOLECULE.atoms_count; j++) {
			double diff_x = kd->charges[atoms_processed + j] - average_calculated_charge;
			double diff_y = MOLECULE.atoms[j].reference_charge - MOLECULE.average_charge;

			cov_xy += diff_x * diff_y;
			cov_xx += diff_x * diff_x;
			cov_yy += diff_y * diff_y;
		}

		kd->per_molecule_stats[i].R = (float) ((cov_xy * cov_xy) / (cov_xx * cov_yy));

		/* Avoid division by zero */
		if(fabs(cov_xx * cov_yy) <= 0.0f)
			bad_molecules++;
		else
			R_sum_molecules += (cov_xy * cov_xy) / (cov_xx * cov_yy);

		atoms_processed += MOLECULE.atoms_count;
		#undef MOLECULE
	}

	kd->full_stats.R = (float) (R_sum_molecules / (ts.molecules_count - bad_molecules));
}


/* Set RMSD for the whole set */
static void set_total_RMSD(struct kappa_data * const kd) {

	assert(kd != NULL);

	int atoms_processed = 0;

	double RMSD_sum_molecules = 0.0;
	for(int i = 0; i < ts.molecules_count; i++) {
		#define MOLECULE ts.molecules[i]

		double diff2_sum_molecule = 0.0;
		for(int j = 0; j < MOLECULE.atoms_count; j++) {
			double diff = kd->charges[atoms_processed + j] - MOLECULE.atoms[j].reference_charge;

			diff2_sum_molecule += fabs(diff) * fabs(diff);
		}

		kd->per_molecule_stats[i].RMSD = (float) sqrt(diff2_sum_molecule / MOLECULE.atoms_count);
		RMSD_sum_molecules += sqrt(diff2_sum_molecule / MOLECULE.atoms_count);

		atoms_processed += MOLECULE.atoms_count;
		#undef MOLECULE
	}

	kd->full_stats.RMSD = (float) (RMSD_sum_molecules / ts.molecules_count);
}


/* Set the average absolute difference for the whole set */
static void set_total_D_avg(struct kappa_data * const kd) {

	assert(kd != NULL);

	int atoms_processed = 0;

	double D_avg_sum_molecules = 0.0;
	for(int i = 0; i < ts.molecules_count; i++) {
		#define MOLECULE ts.molecules[i]

		double D_sum_molecule = 0.0;
		for(int j = 0; j < MOLECULE.atoms_count; j++) {
			double diff = kd->charges[atoms_processed + j] - MOLECULE.atoms[j].reference_charge;

			D_sum_molecule += fabs(diff);
		}

		kd->per_molecule_stats[i].D_avg = (float) (D_sum_molecule / MOLECULE.atoms_count);
		D_avg_sum_molecules += D_sum_molecule / MOLECULE.atoms_count;
		atoms_processed += MOLECULE.atoms_count;
		#undef MOLECULE
	}

	kd->full_stats.D_avg = (float) (D_avg_sum_molecules / ts.molecules_count);
}


/* Set the maximum absolute difference for the whole set */
static void set_total_D_max(struct kappa_data * const kd) {

	assert(kd != NULL);

	int atoms_processed = 0;

	double D_max_sum_molecules = 0.0;
	for(int i = 0; i < ts.molecules_count; i++) {
		#define MOLECULE ts.molecules[i]

		double max_diff_per_molecule = 0.0;
		for(int j = 0; j < MOLECULE.atoms_count; j++) {
			double diff = kd->charges[atoms_processed + j] - MOLECULE.atoms[j].reference_charge;

			if(fabs(diff) > max_diff_per_molecule)
				max_diff_per_molecule = fabs(diff);
		}

		D_max_sum_molecules += max_diff_per_molecule;
		kd->per_molecule_stats[i].D_max = (float) max_diff_per_molecule;

		atoms_processed += MOLECULE.atoms_count;
		#undef MOLECULE
	}

	kd->full_stats.D_max = (float) (D_max_sum_molecules / ts.molecules_count);
}


/* Set Pearson weighted correlation coeff. for the whole set */
static void set_total_R_weighted(struct subset * const ss, struct kappa_data * const kd) {

	int atoms_processed = 0;

	int bad_molecules = 0;
	double Rw_sum_molecules = 0.0;

	for(int i = 0; i < ts.molecules_count; i++) {
		#define MOLECULE ts.molecules[i]

		double avg_wgh_qm_chg = 0.0;
		double avg_wgh_eem_chg = 0.0;
		double weights_sum = 0.0;
		for(int j = 0; j < MOLECULE.atoms_count; j++) {
			double weight = ss->weights[get_atom_type_idx(&MOLECULE.atoms[j])];
			avg_wgh_eem_chg +=  weight * kd->charges[atoms_processed + j];
			avg_wgh_qm_chg += weight * MOLECULE.atoms[j].reference_charge;

			weights_sum += weight;
		}

		avg_wgh_qm_chg /= weights_sum;
		avg_wgh_eem_chg /= weights_sum;

		double cov_xy = 0.0;
		double cov_xx = 0.0;
		double cov_yy = 0.0;

		weights_sum = 0.0;
		for(int j = 0; j < MOLECULE.atoms_count; j++) {
			double weight = ss->weights[get_atom_type_idx(&MOLECULE.atoms[j])];
			double diff_x = kd->charges[atoms_processed + j] - avg_wgh_eem_chg;
			double diff_y = MOLECULE.atoms[j].reference_charge - avg_wgh_qm_chg;
			weights_sum += weight;

			cov_xy += weight * diff_x * diff_y;
			cov_xx += weight * diff_x * diff_x;
			cov_yy += weight * diff_y * diff_y;
		}

		kd->per_molecule_stats[i].R_weighted = (float) ((cov_xy * cov_xy) / (cov_xx * cov_yy));

		/* Avoid division by zero */
		if(fabs(cov_xx * cov_yy) <= 0.0f)
			bad_molecules++;
		else
			Rw_sum_molecules += (cov_xy * cov_xy) / (cov_xx * cov_yy);

		atoms_processed += MOLECULE.atoms_count;
		#undef MOLECULE
	}

	kd->full_stats.R_weighted = (float) (Rw_sum_molecules / (ts.molecules_count - bad_molecules));
}


/* Set Pearson correlation coeff. for each atom type */
static void set_per_at_R(struct kappa_data * const kd) {

	assert(kd != NULL);

	/* Compute starting indices for storing the charges of each molecule. These are
	 * needed to access individual charges. */
	int starts[ts.molecules_count];
	starts[0] = 0;
	for(int i = 1; i < ts.molecules_count; i++)
		starts[i] = starts[i - 1] + ts.molecules[i - 1].atoms_count;

	for(int i = 0; i < ts.atom_types_count; i++) {
		#define AT ts.atom_types[i]

		double avg_eem_chg_per_at = 0.0;
		double avg_qm_chg_per_at = 0.0;

		for(int j = 0; j < ts.atom_types[i].atoms_count; j++) {
			const int molecule_idx = AT.atoms_molecule_idx[j];
			const int atom_idx = AT.atoms_atom_idx[j];

			avg_qm_chg_per_at += ts.molecules[molecule_idx].atoms[atom_idx].reference_charge;
			avg_eem_chg_per_at += kd->charges[starts[molecule_idx] + atom_idx];
		}

		avg_qm_chg_per_at /= AT.atoms_count;
		avg_eem_chg_per_at /= AT.atoms_count;

		double cov_xx = 0.0;
		double cov_xy = 0.0;
		double cov_yy = 0.0;

		for(int j = 0; j < ts.atom_types[i].atoms_count; j++) {
			const int molecule_idx = AT.atoms_molecule_idx[j];
			const int atom_idx = AT.atoms_atom_idx[j];

			double diff_x = kd->charges[starts[molecule_idx] + atom_idx] - avg_eem_chg_per_at;
			double diff_y = ts.molecules[molecule_idx].atoms[atom_idx].reference_charge - avg_qm_chg_per_at;

			cov_xy += diff_x * diff_y;
			cov_xx += diff_x * diff_x;
			cov_yy += diff_y * diff_y;
		}

		kd->per_at_stats[i].R = (float) ((cov_xy * cov_xy) / (cov_xx * cov_yy));
		#undef AT
	}

}


/* Set RMSD for each atom type */
static void set_per_at_RMSD(struct kappa_data * const kd) {

	assert(kd != NULL);

	/* Compute starting indices for storing the charges of each molecule. These are
	 * needed to access individual charges. */
	int starts[ts.molecules_count];
	starts[0] = 0;
	for(int i = 1; i < ts.molecules_count; i++)
		starts[i] = starts[i - 1] + ts.molecules[i - 1].atoms_count;

	for(int i = 0; i < ts.atom_types_count; i++) {
		#define AT ts.atom_types[i]

		double diff2_sum = 0.0;
		for(int j  = 0; j < AT.atoms_count; j++) {
			const int molecule_idx = AT.atoms_molecule_idx[j];
			const int atom_idx = AT.atoms_atom_idx[j];

			double diff = ts.molecules[molecule_idx].atoms[atom_idx].reference_charge - kd->charges[starts[molecule_idx] + atom_idx];
			diff2_sum += diff * diff;
		}

		kd->per_at_stats[i].RMSD = (float) sqrt(diff2_sum / AT.atoms_count);

		#undef AT
	}
}


/* Set average absolute difference for each atom type */
static void set_per_at_D_avg(struct kappa_data * const kd) {

	assert(kd != NULL);

	int atoms_processed = 0;

	double *D_sum_atom_type = (double *) calloc(ts.atom_types_count, sizeof(double));
	if(!D_sum_atom_type)
		EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for statistical data.\n");

	for(int i = 0; i < ts.molecules_count; i++) {
		#define MOLECULE ts.molecules[i]

		for(int j = 0; j < MOLECULE.atoms_count; j++) {
			double diff = kd->charges[atoms_processed + j] - MOLECULE.atoms[j].reference_charge;

			const int at_idx = get_atom_type_idx(&MOLECULE.atoms[j]);

			D_sum_atom_type[at_idx] += fabs(diff);
		}

		atoms_processed += MOLECULE.atoms_count;
		#undef MOLECULE
	}

	for(int i = 0; i < ts.atom_types_count; i++)
		kd->per_at_stats[i].D_avg = (float) D_sum_atom_type[i] / ts.atom_types[i].atoms_count;

	free(D_sum_atom_type);
}


/* Set maximal absolute difference for each atom type */
static void set_per_at_D_max(struct kappa_data * const kd) {

	assert(kd != NULL);

	int atoms_processed = 0;

	/* Reset values left by previous iteration of the Brent's method */
	for(int i = 0; i < ts.atom_types_count; i++)
		kd->per_at_stats[i].D_max = 0.0f;

	for(int i = 0; i < ts.molecules_count; i++) {
		#define MOLECULE ts.molecules[i]

		for(int j = 0; j < MOLECULE.atoms_count; j++) {
			double diff = kd->charges[atoms_processed + j] - MOLECULE.atoms[j].reference_charge;

			const int at_idx = get_atom_type_idx(&MOLECULE.atoms[j]);

			if(fabs(diff) > kd->per_at_stats[at_idx].D_max)
				kd->per_at_stats[at_idx].D_max = (float) fabs(diff);
		}

		atoms_processed += MOLECULE.atoms_count;
		#undef MOLECULE
	}
}


/* Calculate all supported statistics at once */
void calculate_statistics(struct subset * const ss, struct kappa_data * const kd) {

	assert(ss != NULL);
	assert(kd != NULL);

	set_total_R(kd);
	set_total_R_weighted(ss, kd);
	set_total_RMSD(kd);
	set_total_D_avg(kd);
	set_total_D_max(kd);

	/* Calculate per atom type statistics */
	set_per_at_R(kd);
	set_per_at_RMSD(kd);
	set_per_at_D_max(kd);
	set_per_at_D_avg(kd);
}

/* Check for abnormal charge differences */
void check_charges(const struct kappa_data * const kd) {

	assert(kd != NULL);

	int bad_molecules = 0;

	for(int i = 0; i < ts.molecules_count; i++) {
		if(kd->per_molecule_stats[i].R < WARN_MIN_R ||
		   kd->per_molecule_stats[i].RMSD > WARN_MAX_RMSD ||
		   kd->per_molecule_stats[i].D_avg > WARN_MAX_D_AVG ||
		   kd->per_molecule_stats[i].D_max > WARN_MAX_D_MAX) {
			fprintf(stderr, "Warning: Abnormal values for molecule %s\n", ts.molecules[i].name);
			fprintf(stderr, "R: %6.4f  RMSD: %4.2e  D_avg: %4.2e  D_max: %4.2e\n\n",
				kd->per_molecule_stats[i].R, kd->per_molecule_stats[i].RMSD,
				kd->per_molecule_stats[i].D_avg, kd->per_molecule_stats[i].D_max);

			bad_molecules++;
		}

	}
	fprintf(stderr, "Check charges: %d molecules with abnormal statistics found.\n", bad_molecules);
}
