/*
 * NEEMP - structures.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "neemp.h"
#include "settings.h"
#include "structures.h"

extern const struct settings s;
extern struct training_set ts;

static void a_destroy(struct atom * const a);
static double rdist(const struct atom * const a1, const struct atom * const a2);
static void m_calculate_rdists(struct molecule * const m);
static void m_calculate_avg_electronegativity(struct molecule * const m);
static void m_calculate_charge_stats(struct molecule * const m);
static void fill_atom_types(void);
static void discard_molecules_without_charges(void);
static void calculate_y(void);

/* Symbols for chemical elements */
static const char * const elems[] = {"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr"};

/* Electronegativities of the chemical elements */
static const float electronegativies[] = {2.2, 0, 0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98, 0, 0.93, 1.31, 1.61, 1.9, 2.19, 2.58, 3.16, 0, 0.82, 1, 1.36, 1.54, 1.63, 1.66, 1.55, 1.83, 1.88, 1.91, 1.9, 1.65, 1.81, 2.01, 2.18, 2.55, 2.96, 3, 0.82, 0.95, 1.22, 1.33, 1.6, 2.16, 1.9, 2.2, 2.28, 2.2, 1.93, 1.69, 1.78, 1.96, 2.05, 2.1, 2.66, 2.6, 0.79, 0.89, 1.1, 1.12, 1.13, 1.14, 1.13, 1.17, 1.2, 1.2, 1.2, 1.22, 1.23, 1.24, 1.25, 1.1, 1.27, 1.3, 1.5, 2.36, 1.9, 2.2, 2.2, 2.28, 2.54, 2, 2.04, 2.33, 2.02, 2, 2.2, 0, 0.7, 0.9, 1.1, 1.3, 1.5, 1.38, 1.36, 1.28, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3};

/* Convert atomic symbol to atomic number */
int convert_symbol_to_Z(const char * const symb) {

	assert(symb != NULL);

	for(int i = 0; i < 100; i++)
		if(!strcmp(symb, elems[i]))
			return i + 1;

	/* Not found */
	return 0;
}

/* Convert atomic number to symbol */
const char *convert_Z_to_symbol(int Z) {

	if(0 < Z && Z < 100)
		return elems[Z - 1];
	else
		return "??";
}

/* Destroy contents of an atom */
static void a_destroy(struct atom * const a) {

	assert(a != NULL);

	free(a->rdists);
}

/* Calculates reciprocal distance between two atoms */
static double rdist(const struct atom * const a1, const struct atom * const a2) {

	assert(a1 != NULL);
	assert(a2 != NULL);

	double dx = a1->position[0] - a2->position[0];
	double dy = a1->position[1] - a2->position[1];
	double dz = a1->position[2] - a2->position[2];

	return 1.0 / sqrt(dx * dx + dy * dy + dz * dz);
}

/* Calculate reciprocal distances for all atoms in the molecule */
static void m_calculate_rdists(struct molecule * const m) {

	assert(m != NULL);

	for(int i = 0; i < m->atoms_count; i++) {
		m->atoms[i].rdists[i] = 0.0;
		for(int j = i + 1; j < m->atoms_count; j++) {
			double rd = rdist(&m->atoms[i], &m->atoms[j]);
			/* Distance is symmetrical */
			m->atoms[i].rdists[j] = rd;
			m->atoms[j].rdists[i] = rd;
		}
	}
}

/* Destroy content of the molecule */
void m_destroy(struct molecule * const m) {

	assert(m != NULL);

	for(int i = 0; i < m->atoms_count; i++)
		a_destroy(&m->atoms[i]);

	free(m->atoms);
	free(m->name);
}

/* Destroy content of the atom type */
void at_destroy(struct atom_type * const at) {

	assert(at != NULL);

	free(at->atoms_molecule_idx);
	free(at->atoms_atom_idx);
}

/* Destroy content of the training set */
void ts_destroy(void) {

	for(int i = 0; i < ts.molecules_count; i++)
		m_destroy(&ts.molecules[i]);

	for(int i = 0; i < ts.atom_types_count; i++)
		at_destroy(&ts.atom_types[i]);

	free(ts.molecules);
	free(ts.atom_types);
}

/* Remove molecules for which the charges haven't been loaded */
static void discard_molecules_without_charges(void) {

	int i = 0;
	while(i < ts.molecules_count) {
		if(!ts.molecules[i].charges_loaded) {
			/* Destroy molecule without charges; fill its space with the last one */
			ts.atoms_count -= ts.molecules[i].atoms_count;
			m_destroy(&ts.molecules[i]);
			ts.molecules_count--;
			/* Check if the last molecule was discarded */
			if(i != ts.molecules_count)
				ts.molecules[i] = ts.molecules[ts.molecules_count];
		} else
			i++;

	}

	/* Shrink molecules array */
	ts.molecules = (struct molecule *) realloc(ts.molecules, sizeof(struct molecule) * ts.molecules_count);
}

/* Fill atom type structures */
static void fill_atom_types(void) {

	ts.atom_types = (struct atom_type *) malloc(sizeof(struct atom_type) * MAX_ATOM_TYPES);
	ts.atom_types_count = 0;

	/* Find number of atom types present in the training set with appropriate atom counts */
	for(int i = 0; i < ts.molecules_count; i++)
		for(int j = 0; j < ts.molecules[i].atoms_count; j++) {
			#define ATOM ts.molecules[i].atoms[j]
			int found = 0;
			for(int k = 0; k < ts.atom_types_count; k++)
				if(ATOM.Z == ts.atom_types[k].Z && ATOM.bond_order == ts.atom_types[k].bond_order) {
					found = 1;
					ts.atom_types[k].atoms_count++;
					break;
				}

			if(!found) {
				/* Create new atom type */
				ts.atom_types[ts.atom_types_count].Z = ATOM.Z;
				ts.atom_types[ts.atom_types_count].bond_order = ATOM.bond_order;
				ts.atom_types[ts.atom_types_count].atoms_count = 1;
				ts.atom_types_count++;
			}
			#undef ATOM
		}

	/* Shrink atom types array */
	ts.atom_types = (struct atom_type *) realloc(ts.atom_types, sizeof(struct atom_type) * ts.atom_types_count);

	/* Allocate space for atom type indices */
	for(int i = 0; i < ts.atom_types_count; i++) {
		ts.atom_types[i].atoms_molecule_idx = (int *) malloc(sizeof(int) * ts.atom_types[i].atoms_count);
		ts.atom_types[i].atoms_atom_idx = (int *) malloc(sizeof(int) * ts.atom_types[i].atoms_count);
		if(!ts.atom_types[i].atoms_molecule_idx || !ts.atom_types[i].atoms_atom_idx)
			EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for atom types.\n");
	}

	/* Reset atom counts */
	for(int i = 0; i < ts.atom_types_count; i++)
		ts.atom_types[i].atoms_count = 0;

	/* Fill atom types with indices pointing to atoms of that kind */
	for(int i = 0; i < ts.molecules_count; i++)
		for(int j = 0; j < ts.molecules[i].atoms_count; j++) {
			#define ATOM ts.molecules[i].atoms[j]
			for(int k = 0; k < ts.atom_types_count; k++) {
				#define AT ts.atom_types[k]
				/* Note that for every atom there exists k which satisfies following condition */
				if(ATOM.Z == AT.Z && ATOM.bond_order == AT.bond_order) {
					AT.atoms_molecule_idx[AT.atoms_count] = i;
					AT.atoms_atom_idx[AT.atoms_count] = j;
					AT.atoms_count++;
					break;
				}
				#undef AT
			}
			#undef ATOM
		}
}

/* Calculate average electronegativity of a molecule (harmonic mean) */
static void m_calculate_avg_electronegativity(struct molecule * const m) {

	assert(m != NULL);

	double hsum = 0.0;
	for(int i = 0; i < m->atoms_count; i++)
		hsum += 1.0 / electronegativies[m->atoms[i].Z - 1];

	m->electronegativity = (float) (m->atoms_count / hsum);
}

/* Calculate sum and average charge of atoms in the molecule */
static void m_calculate_charge_stats(struct molecule * const m) {

	assert(m != NULL);

	double sum = 0.0;
	for(int i = 0; i < m->atoms_count; i++)
		sum += m->atoms[i].reference_charge;

	m->sum_of_charges = (float) sum;
	m->average_charge = m->sum_of_charges / m->atoms_count;
}

/* Calculate y sums for all molecules */
static void calculate_y(void) {

	for(int i = 0; i < ts.molecules_count; i++)
		for(int j = 0; j < ts.molecules[i].atoms_count; j++) {
			#define ATOM(x) ts.molecules[i].atoms[x]
			ATOM(j).y = 0.0;
			for(int k = 0; k < ts.molecules[i].atoms_count; k++) {
				if(j == k)
					continue;
				ATOM(j).y = ATOM(k).reference_charge * ATOM(j).rdists[k];
			}
			#undef ATOM
		}
}

/* Find atom type for a particular atom */
int get_atom_type_idx(const struct atom * const a) {

	assert(a != NULL);

	for(int i = 0; i < ts.atom_types_count; i++)
		if(a->Z == ts.atom_types[i].Z && a->bond_order == ts.atom_types[i].bond_order)
			return i;

	/* We should not get here! */
	assert(0);
	return NOT_FOUND;
}

/* Do some preprocessing to simplify things later on */
void preprocess_molecules(void) {

	/* Remove molecule for which the charges are not present */
	if(s.mode == MODE_PARAMS)
		discard_molecules_without_charges();

	/* Calculate reciprocal distances of atoms for all molecules */
	for(int i = 0; i < ts.molecules_count; i++)
		m_calculate_rdists(&ts.molecules[i]);

	/* Calculate average electronegativies */
	for(int i = 0; i < ts.molecules_count; i++)
		m_calculate_avg_electronegativity(&ts.molecules[i]);

	/* Calculate sum and average of the charges in the molecule */
	if(s.mode == MODE_PARAMS)
		for(int i = 0; i < ts.molecules_count; i++)
			m_calculate_charge_stats(&ts.molecules[i]);

	/* Calculate auxiliary sum */
	if(s.mode == MODE_PARAMS)
		calculate_y();

	/* Finally, fill indices to atoms of a particular kind */
	fill_atom_types();

}

/* Prints information about the training set */
void ts_info(void) {

	printf("Training set info\n");

	printf("Molecules: %5d  Atoms: %8d  Atom types: %2d\n", ts.molecules_count, ts.atoms_count, ts.atom_types_count);
	printf("Atom type    # atoms      %% atoms \n");

	for(int i = 0; i < ts.atom_types_count; i++) {
		#define AT ts.atom_types[i]
		printf(" %2s(%1d)     %8d      %6.3f\n", convert_Z_to_symbol(AT.Z), AT.bond_order, AT.atoms_count, 100.0f * (float) AT.atoms_count / ts.atoms_count);
		#undef AT
	}
}
