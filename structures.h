/*
 * NEEMP - structures.h
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#ifndef __STRUCTURES_H__
#define __STRUCTURES_H__

struct atom {

	int Z;				/* atomic number */
	int bond_order;
	float position[3];		/* x, y, z position */
	float reference_charge;
};

int convert_symbol_to_Z(const char * const symb);
const char *convert_Z_to_symbol(int Z);

struct molecule {

	char *name;
	int atoms_count;
	struct atom *atoms;

	/* Auxiliary variables */

	int charges_loaded;
	float electronegativity;
	float average_charge;
	float sum_of_charges;
};

void m_destroy(struct molecule * const m);

struct atom_type {

	int Z;
	int bond_order;

	int atoms_count;

	/* the pair (atoms_molecule_idx[i], atoms_atom_idx[i]) uniquely
	 * identifies particular atom i of this atom type */
	int *atoms_molecule_idx;
	int *atoms_atom_idx;
};

void at_destroy(struct atom_type * const at);

/* Note that there is only one instance of training_set */
struct training_set {

	struct molecule *molecules;
	struct atom_type *atom_types;

	int atoms_count;
	int atom_types_count;
	int molecules_count;
};

void ts_destroy(void);

void preprocess_molecules(void);

#endif /* __STRUCTURES_H__ */
