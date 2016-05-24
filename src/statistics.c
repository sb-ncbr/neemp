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


static int compare(const void *p1, const void *p2);
static void adjust_ranks_via_pointers(float **array, int n);

static void set_total_Spearman(struct kappa_data * const kd);
static void set_total_R(struct kappa_data * const kd);
static void set_total_R_w(struct kappa_data *const kd);
static void set_total_RMSD(struct kappa_data * const kd);
static void set_total_RMSD_avg(struct kappa_data * const kd);
static void set_total_D_avg(struct kappa_data * const kd);
static void set_total_D_max(struct kappa_data * const kd);

static void set_per_at_R_R2(struct kappa_data * const kd);
static void set_per_at_RMSD(struct kappa_data * const kd);
static void set_per_at_D_avg(struct kappa_data * const kd);
static void set_per_at_D_max(struct kappa_data * const kd);


/* Compare two floats via pointers */
static int compare(const void *p1, const void *p2) {

	float a = **(float **) p1;
	float b = **(float **) p2;

	if(a > b)
		return 1;
	if(a < b)
		return -1;

	return 0;
}


/* Set ranks for Spearman correlation coeff */
static void adjust_ranks_via_pointers(float **array, int n) {

	int latest_rank = 1;
	int i = 0;
	while(i < n) {
		int j = 1;
		while(i + j < n && fabsf(*array[i] - *array[i + j]) < 0.00001f)
			j++;

		float rank = (float) (2.0f * latest_rank + j - 1) / 2.0f;
		for(int k = 0; k < j; k++)
			*array[i + k] = rank;

		latest_rank += j;
		i += j;
	}
}


/* Set total weighted correlation computed of individual Pearson's coeff per atom type */
static void set_total_R_w(struct kappa_data * const kd) {

	assert(kd != NULL);

	double weighted_corr_sum = 0.0;

	for(int i = 0; i < ts.atom_types_count; i++) {
			double weight = pow(s.rw, kd->per_at_stats[i].R);
			weighted_corr_sum += weight * kd->per_at_stats[i].R;
	}

	//consider total R
	weighted_corr_sum += 3 * (pow(s.rw, kd->full_stats.R2)) * kd->full_stats.R2 - kd->full_stats.RMSD/3;

	/* Normalize the results */
	kd->full_stats.R_w = (float) (weighted_corr_sum / (ts.atom_types_count * s.rw + 3 * s.rw));
}


/* Set total Spearman correlation coeff. for the whole set */
static void set_total_Spearman(struct kappa_data * const kd) {

	assert(kd != NULL);

	int atoms_processed = 0;
	int bad_molecules = 0;

	double spearman_sum_molecules = 0.0;

	for(int i = 0; i < ts.molecules_count; i++) {
		#define MOLECULE ts.molecules[i]

		float *reference_data = (float *) calloc(MOLECULE.atoms_count, sizeof(float));
		float *calculated_data = (float *) calloc(MOLECULE.atoms_count, sizeof(float));
		float **reference_data_pointers = (float **) calloc(MOLECULE.atoms_count, sizeof(float *));
		float **calculated_data_pointers = (float **) calloc(MOLECULE.atoms_count, sizeof(float *));
		if(!reference_data || !calculated_data || !calculated_data_pointers || !reference_data_pointers)
			EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for Spearman correlation computation.");

		memcpy(calculated_data, &kd->charges[atoms_processed], sizeof(float) * MOLECULE.atoms_count);
		for(int j = 0; j < MOLECULE.atoms_count; j++)
			reference_data[j] = MOLECULE.atoms[j].reference_charge;

		/* Set pointers to the data */
		for(int j = 0; j < MOLECULE.atoms_count; j++) {
			reference_data_pointers[j] = &reference_data[j];
			calculated_data_pointers[j] = &calculated_data[j];
		}

		qsort(reference_data_pointers, MOLECULE.atoms_count, sizeof(float *), compare);
		qsort(calculated_data_pointers, MOLECULE.atoms_count, sizeof(float *), compare);

		adjust_ranks_via_pointers(reference_data_pointers, MOLECULE.atoms_count);
		adjust_ranks_via_pointers(calculated_data_pointers, MOLECULE.atoms_count);

		/* Use Pearson correlation between computed ranks */
		double average_calculated_rank = 0.0;
		double average_reference_rank = 0.0;
		for(int j = 0; j < MOLECULE.atoms_count; j++) {
			average_reference_rank += reference_data[j];
			average_calculated_rank += calculated_data[j];
		}

		average_reference_rank /= MOLECULE.atoms_count;
		average_calculated_rank /= MOLECULE.atoms_count;

		double cov_xy = 0.0;
		double cov_xx = 0.0;
		double cov_yy = 0.0;

		for(int j = 0; j < MOLECULE.atoms_count; j++) {
			double diff_x = calculated_data[j] - average_calculated_rank;
			double diff_y = reference_data[j] - average_reference_rank;

			cov_xy += diff_x * diff_y;
			cov_xx += diff_x * diff_x;
			cov_yy += diff_y * diff_y;
		}

		kd->per_molecule_stats[i].spearman = (float) (cov_xy / sqrt(cov_xx * cov_yy));

		/* Avoid division by zero */
		if(fabs(cov_xx * cov_yy) <= 0.0f)
			bad_molecules++;
		else
			spearman_sum_molecules += cov_xy / sqrt(cov_xx * cov_yy);

		free(reference_data);
		free(calculated_data);
		free(reference_data_pointers);
		free(calculated_data_pointers);

		atoms_processed += MOLECULE.atoms_count;
		#undef MOLECULE
	}

	kd->full_stats.spearman = (float) (spearman_sum_molecules / (ts.molecules_count - bad_molecules));
}

/* Set total Pearson correlation coeff for the whole set */
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

		kd->per_molecule_stats[i].R = (float) (cov_xy / sqrt(cov_xx * cov_yy));

		/* Avoid division by zero */
		if(fabs(cov_xx * cov_yy) <= 0.0f)
			bad_molecules++;
		else
			R_sum_molecules += cov_xy / sqrt(cov_xx * cov_yy);

		atoms_processed += MOLECULE.atoms_count;
		#undef MOLECULE
	}

	kd->full_stats.R = (float) (R_sum_molecules / (ts.molecules_count - bad_molecules));
}

/* Set total Pearson correlation coeff for the whole set */
static void set_total_R2(struct kappa_data * const kd) {

	assert(kd != NULL);

	double R2_sum_molecules = 0.0;

	for(int i = 0; i < ts.molecules_count; i++) {
		double molecule_R = kd->per_molecule_stats[i].R;
		kd->per_molecule_stats[i].R2 = (float) (molecule_R * molecule_R);
		R2_sum_molecules += molecule_R * molecule_R;
	}

	kd->full_stats.R2 = (float) R2_sum_molecules / ts.molecules_count;
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

/* Set RMSD_avg for the whole set */
static void set_total_RMSD_avg(struct kappa_data * const kd) {
    assert(kd!=NULL);

    double RMSD_sum_atom_types = 0.0;
    for (int i=0; i < ts.atom_types_count; i++) {
        RMSD_sum_atom_types += kd->per_at_stats[i].RMSD;
    }
    kd->full_stats.RMSD_avg = (float) (RMSD_sum_atom_types / ts.atom_types_count);
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

/* Set Pearson correlation coeff. for each atom type */
static void set_per_at_R_R2(struct kappa_data * const kd) {

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

		for(int j = 0; j < AT.atoms_count; j++) {
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

		kd->per_at_stats[i].R = (float) ((cov_xy) / sqrt(cov_xx * cov_yy));
		kd->per_at_stats[i].R2 = (float) ((cov_xy * cov_xy) / (cov_xx * cov_yy));
		#undef AT
	}

}


/* Set Spearman correlation coeff for each atom type */
static void set_per_at_Spearman(struct kappa_data * const kd) {

	assert(kd != NULL);

	/* Compute starting indices for storing the charges of each molecule. These are
	 * needed to access individual charges. */
	int starts[ts.molecules_count];
	starts[0] = 0;
	for(int i = 1; i < ts.molecules_count; i++)
		starts[i] = starts[i - 1] + ts.molecules[i - 1].atoms_count;

	for(int i = 0; i < ts.atom_types_count; i++) {
		#define AT ts.atom_types[i]
		float *reference_data = (float *) calloc(AT.atoms_count, sizeof(float));
		float *calculated_data = (float *) calloc(AT.atoms_count, sizeof(float));
		float **calculated_data_pointers = (float **) calloc(AT.atoms_count, sizeof(float *));
		float **reference_data_pointers = (float **) calloc(AT.atoms_count, sizeof(float *));

		if(!reference_data || !calculated_data || !calculated_data_pointers || !reference_data_pointers)
			EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for Spearman correlation computation.");

		for(int j = 0; j < AT.atoms_count; j++) {
			const int molecule_idx = AT.atoms_molecule_idx[j];
			const int atom_idx = AT.atoms_atom_idx[j];

			calculated_data[j] = kd->charges[starts[molecule_idx] + atom_idx];
			reference_data[j] = ts.molecules[molecule_idx].atoms[atom_idx].reference_charge;
		}

		/* Set pointers to the data */
		for(int j = 0; j < AT.atoms_count; j++) {
			reference_data_pointers[j] = &reference_data[j];
			calculated_data_pointers[j] = &calculated_data[j];
		}

		qsort(reference_data_pointers, AT.atoms_count, sizeof(float *), compare);
		qsort(calculated_data_pointers, AT.atoms_count, sizeof(float *), compare);

		adjust_ranks_via_pointers(reference_data_pointers, AT.atoms_count);
		adjust_ranks_via_pointers(calculated_data_pointers, AT.atoms_count);

		/* Use Pearson correlation between computed ranks */
		double average_calculated_rank = 0.0;
		double average_reference_rank = 0.0;
		for(int j = 0; j < AT.atoms_count; j++) {
			average_reference_rank += reference_data[j];
			average_calculated_rank += calculated_data[j];
		}

		average_reference_rank /= AT.atoms_count;
		average_calculated_rank /= AT.atoms_count;

		double cov_xy = 0.0;
		double cov_xx = 0.0;
		double cov_yy = 0.0;

		for(int j = 0; j < AT.atoms_count; j++) {
			double diff_x = calculated_data[j] - average_calculated_rank;
			double diff_y = reference_data[j] - average_reference_rank;

			cov_xy += diff_x * diff_y;
			cov_xx += diff_x * diff_x;
			cov_yy += diff_y * diff_y;
		}

		kd->per_at_stats[i].spearman = (float) (cov_xy / sqrt(cov_xx * cov_yy));

		free(reference_data);
		free(calculated_data);
		free(reference_data_pointers);
		free(calculated_data_pointers);
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

	/* Calculate total statistics */
	set_total_Spearman(kd);
	set_total_R(kd);
	set_total_R2(kd);
	set_total_RMSD(kd);
	set_total_D_avg(kd);
	set_total_D_max(kd);

	/* Calculate per atom type statistics */
	set_per_at_R_R2(kd);
	set_per_at_Spearman(kd);
	set_per_at_RMSD(kd);
	set_per_at_D_max(kd);
	set_per_at_D_avg(kd);

	/* Computed from per atom type stats, needs to go last */
	set_total_R_w(kd);
	set_total_RMSD_avg(kd);
}

/* Calculate statistics according to set sort type */
void calculate_statistics_by_sort_mode(struct kappa_data* kd) {
	assert(kd != NULL);
	switch (s.sort_by) {
		case SORT_R:
		case SORT_R2:
			set_total_R(kd);
			set_total_R2(kd);
			set_per_at_R_R2(kd);
			break;
		case SORT_RMSD:
		case SORT_RMSD_AVG:
			set_total_RMSD(kd);
			set_per_at_RMSD(kd);
			set_total_RMSD_avg(kd);
			break;
		case SORT_SPEARMAN:
			set_total_Spearman(kd);
			set_per_at_Spearman(kd);
			break;
		case SORT_D_AVG:
			set_total_D_avg(kd);
			set_per_at_D_avg(kd);
			break;
		case SORT_D_MAX:
			set_total_D_max(kd);
			set_per_at_D_max(kd);
			break;
		case SORT_RW:
			set_total_R(kd);
			set_per_at_R_R2(kd);
			set_total_RMSD(kd);
			set_per_at_RMSD(kd);
			set_total_R_w(kd);
			break;
		default:
			break;

	}
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
