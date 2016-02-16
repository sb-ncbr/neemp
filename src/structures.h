/*
 * NEEMP - structures.h
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013, 2014
 *
 * */

#ifndef __STRUCTURES_H__
#define __STRUCTURES_H__

struct atom {

	int Z;				/* atomic number */
	int bond_order;
	float position[3];		/* x, y, z position */
	float reference_charge;

	double *rdists;
	double y;			/* sum of charges/distance over all atoms in the molecule */
};

int convert_symbol_to_Z(const char * const symb);
const char *convert_Z_to_symbol(int Z);
double rdist(const struct atom * const a1, const struct atom * const a2);

struct molecule {

	char *name;
	int atoms_count;
	struct atom *atoms;

	/* Auxiliary variables */
	int is_valid;
	int has_charges;
	int has_parameters;
	float electronegativity;
	float average_charge;
	float sum_of_charges;
};

void m_destroy(struct molecule * const m);
void get_sum_formula(const struct molecule * const m, char * const buff, int n);

struct atom_type {

	int Z;
	int bond_order;

	int atoms_count;
	int molecules_count;

	/* the pair (atoms_molecule_idx[i], atoms_atom_idx[i]) uniquely
	 * identifies particular atom i of this atom type */
	int *atoms_molecule_idx;
	int *atoms_atom_idx;

	int has_parameters;
};

int get_atom_type_idx(const struct atom * const a);
int get_atom_type_idx_from_text(const char * const str);
void at_destroy(struct atom_type * const at);
void at_format_text(const struct atom_type * const at, char * const buff);

/* Note that there is only one instance of training_set */
struct training_set {

	struct molecule *molecules;
	struct atom_type *atom_types;

	int atoms_count;
	int atom_types_count;
	int molecules_count;
};

/* Electron affinity if the chemical elements */
static const float affinities[] = {0.0, 0.75420375, 0.0, 0.618049, 0.0, 0.279723, 1.262118, -0.07, 1.461112, 3.4011887, 0.0, 0.547926, 0.0, 0.43283, 1.389521, 0.7465, 2.0771029, 3.612724, 0.0, 0.501459, 0.02455, 0.188, 0.084, 0.525, 0.67584, 0.0, 0.151, 0.6633, 1.15716, 1.23578, 0.0, 0.41, 1.232712, 0.814, 2.02067, 3.363588, 0.0, 0.485916, 0.05206, 0.307, 0.426, 0.893, 0.7472, 0.55, 1.04638, 1.14289, 0.56214, 1.30447, 0.0, 0.404, 1.112066, 1.047401, 1.970875, 3.059038, 0.0, 0.471626, 0.14462, 0.47, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, 0.322, 0.815, 0.15, 1.0778, 1.56436, 2.1251, 2.30861, 0.0, 0.377, 0.364, 0.942363, 1.9, 2.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

/* Ion energy of the chemical atoms */
static const float ionenergies[] = {0.0, 13.5984, 24.5874, 5.3917, 9.3227, 8.298, 11.2603, 14.5341, 13.6181, 17.4228, 21.5645, 5.1391, 7.6462, 5.9858, 8.1517, 10.4867, 10.36, 12.9676, 15.7596, 4.3407, 6.1132, 6.5615, 6.8281, 6.7462, 6.7665, 7.434, 7.9024, 7.881, 7.6398, 7.7264, 9.3942, 5.9993, 7.8994, 9.7886, 9.7524, 11.8138, 13.9996, 4.1771, 5.6949, 6.2173, 6.6339, 6.7589, 7.0924, 7.28, 7.3605, 7.4589, 8.3369, 7.5762, 8.9938, 5.7864, 7.3439, 8.6084, 9.0096, 10.4513, 12.1298, 3.8939, 5.2117, 5.5769, 5.5387, 5.473, 5.525, 5.582, 5.6437, 5.6704, 6.1498, 5.8638, 5.9389, 6.0215, 6.1077, 6.1843, 6.2542, 5.4259, 6.8251, 7.5496, 7.864, 7.8335, 8.4382, 8.967, 8.9588, 9.2255, 10.4375, 6.1082, 7.4167, 7.2855, 8.414, 0, 10.7485, 4.0727, 5.2784, 5.17, 6.3067, 5.89, 6.1941, 6.2657, 6.026, 5.9738, 5.9914, 6.1979, 6.2817, 6.42, 6.5, 6.58, 6.65, 4.9, 6.0, 0, 0, 0, 0, 0};
void ts_destroy(void);
void ts_info(void);

void preprocess_molecules(void);
void discard_invalid_molecules_or_without_charges_or_parameters(void);

#endif /* __STRUCTURES_H__ */
