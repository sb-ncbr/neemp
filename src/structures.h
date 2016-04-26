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

	char type_string[10];

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
	int has_atom_types;
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

	char type_string[10];

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

void ts_destroy(void);
void ts_info(void);

void preprocess_molecules(void);
void discard_invalid_molecules_or_without_charges_or_parameters(void);

#endif /* __STRUCTURES_H__ */
