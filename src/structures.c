/*
 * NEEMP - structures.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013, 2014
 *
 * */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "config.h"
#include "neemp.h"
#include "settings.h"
#include "structures.h"

extern const struct settings s;
extern struct training_set ts;

static void a_destroy(struct atom * const a);
static void m_calculate_rdists(struct molecule * const m);
static void m_calculate_avg_electronegativity(struct molecule * const m);
static void m_calculate_charge_stats(struct molecule * const m);
static void fill_atom_types(void);
static void list_molecules_without_charges(void);
static void list_molecules_without_parameters(void);
static void list_invalid_molecules(void);
static void calculate_y(void);
static void at_fill_from_atom(struct atom_type * const at, const struct atom * const a);
static int at_compare_against_atom(const struct atom_type * const at, const struct atom * const a);

/* Symbols for chemical elements */
static const char * const elems[] = {"??", "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr"};

/* Electronegativities of the chemical elements */
static const float electronegativies[] = {0.0, 2.2, 0, 0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98, 0, 0.93, 1.31, 1.61, 1.9, 2.19, 2.58, 3.16, 0, 0.82, 1, 1.36, 1.54, 1.63, 1.66, 1.55, 1.83, 1.88, 1.91, 1.9, 1.65, 1.81, 2.01, 2.18, 2.55, 2.96, 3, 0.82, 0.95, 1.22, 1.33, 1.6, 2.16, 1.9, 2.2, 2.28, 2.2, 1.93, 1.69, 1.78, 1.96, 2.05, 2.1, 2.66, 2.6, 0.79, 0.89, 1.1, 1.12, 1.13, 1.14, 1.13, 1.17, 1.2, 1.2, 1.2, 1.22, 1.23, 1.24, 1.25, 1.1, 1.27, 1.3, 1.5, 2.36, 1.9, 2.2, 2.2, 2.28, 2.54, 2, 2.04, 2.33, 2.02, 2, 2.2, 0, 0.7, 0.9, 1.1, 1.3, 1.5, 1.38, 1.36, 1.28, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3};

/* Electron affinity if the chemical elements */
const float affinities[] = {0.0, 0.75420375, 0.0, 0.618049, 0.0, 0.279723, 1.262118, -0.07, 1.461112, 3.4011887, 0.0, 0.547926, 0.0, 0.43283, 1.389521, 0.7465, 2.0771029, 3.612724, 0.0, 0.501459, 0.02455, 0.188, 0.084, 0.525, 0.67584, 0.0, 0.151, 0.6633, 1.15716, 1.23578, 0.0, 0.41, 1.232712, 0.814, 2.02067, 3.363588, 0.0, 0.485916, 0.05206, 0.307, 0.426, 0.893, 0.7472, 0.55, 1.04638, 1.14289, 0.56214, 1.30447, 0.0, 0.404, 1.112066, 1.047401, 1.970875, 3.059038, 0.0, 0.471626, 0.14462, 0.47, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, 0.322, 0.815, 0.15, 1.0778, 1.56436, 2.1251, 2.30861, 0.0, 0.377, 0.364, 0.942363, 1.9, 2.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

/* Ion energy of the chemical atoms */
const float ionenergies[] = {0.0, 13.5984, 24.5874, 5.3917, 9.3227, 8.298, 11.2603, 14.5341, 13.6181, 17.4228, 21.5645, 5.1391, 7.6462, 5.9858, 8.1517, 10.4867, 10.36, 12.9676, 15.7596, 4.3407, 6.1132, 6.5615, 6.8281, 6.7462, 6.7665, 7.434, 7.9024, 7.881, 7.6398, 7.7264, 9.3942, 5.9993, 7.8994, 9.7886, 9.7524, 11.8138, 13.9996, 4.1771, 5.6949, 6.2173, 6.6339, 6.7589, 7.0924, 7.28, 7.3605, 7.4589, 8.3369, 7.5762, 8.9938, 5.7864, 7.3439, 8.6084, 9.0096, 10.4513, 12.1298, 3.8939, 5.2117, 5.5769, 5.5387, 5.473, 5.525, 5.582, 5.6437, 5.6704, 6.1498, 5.8638, 5.9389, 6.0215, 6.1077, 6.1843, 6.2542, 5.4259, 6.8251, 7.5496, 7.864, 7.8335, 8.4382, 8.967, 8.9588, 9.2255, 10.4375, 6.1082, 7.4167, 7.2855, 8.414, 0, 10.7485, 4.0727, 5.2784, 5.17, 6.3067, 5.89, 6.1941, 6.2657, 6.026, 5.9738, 5.9914, 6.1979, 6.2817, 6.42, 6.5, 6.58, 6.65, 4.9, 6.0, 0, 0, 0, 0, 0};


/* Convert atomic symbol to atomic number */
int convert_symbol_to_Z(const char * const symb) {

	assert(symb != NULL);

	for(int i = 1; i < 103; i++)
		if(!strcasecmp(symb, elems[i]))
			return i;

	/* Check for deuterium and tritium */
	if(!strcmp(symb, "D") || !strcmp(symb, "T"))
		return 1;

	/* Check for unspecified atom by Mol file format */
	if(!strcmp(symb, "A") || !strcmp(symb, "*") || !strcmp(symb, "R"))
		return 0;

	/* Not found */
	return 0;
}

/* Convert atomic number to symbol */
const char *convert_Z_to_symbol(int Z) {

	if(0 < Z && Z < 103)
		return elems[Z];
	else
		return "??";
}

/* Destroy contents of an atom */
static void a_destroy(struct atom * const a) {

	assert(a != NULL);

	if(s.mode == MODE_PARAMS)
		free(a->rdists);
}

/* Calculates reciprocal distance between two atoms */
double rdist(const struct atom * const a1, const struct atom * const a2) {

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
				if(at_compare_against_atom(&ts.atom_types[k], &ATOM)) {
					found = 1;
					ts.atom_types[k].atoms_count++;
					break;
				}

			if(!found) {
				/* Create new atom type */
				at_fill_from_atom(&ts.atom_types[ts.atom_types_count], &ATOM);
				ts.atom_types[ts.atom_types_count].atoms_count = 1;
				ts.atom_types[ts.atom_types_count].has_parameters = 0;
				ts.atom_types_count++;
			}
			#undef ATOM
		}

	if(ts.atom_types_count > MAX_ATOM_TYPES)
		EXIT_ERROR(RUN_ERROR, "Maximum number of atom types (%d) reached. "
				      "Increase value of MAX_ATOM_TYPES in config.h and recompile NEEMP.\n",
				      MAX_ATOM_TYPES);

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
				if(at_compare_against_atom(&AT, &ATOM)) {
					AT.atoms_molecule_idx[AT.atoms_count] = i;
					AT.atoms_atom_idx[AT.atoms_count] = j;
					AT.atoms_count++;
					break;
				}
				#undef AT
			}
			#undef ATOM
		}

	/* Calculate number of molecules which contain particular atom type */
	for(int i = 0; i < ts.atom_types_count; i++) {
		#define AT ts.atom_types[i]
		int *molecule_indices = (int *) calloc(ts.molecules_count, sizeof(int));
		AT.molecules_count = 0;
		for(int j = 0; j < AT.atoms_count; j++)
			molecule_indices[AT.atoms_molecule_idx[j]] = 1;

		for(int j = 0; j < ts.molecules_count; j++)
			if(molecule_indices[j])
				AT.molecules_count++;
		#undef AT
		free(molecule_indices);
	}
}

/* Calculate average electronegativity of a molecule (harmonic mean) */
static void m_calculate_avg_electronegativity(struct molecule * const m) {

	assert(m != NULL);

	double hsum = 0.0;
	for(int i = 0; i < m->atoms_count; i++)
		hsum += 1.0 / electronegativies[m->atoms[i].Z];

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
				ATOM(j).y += ATOM(k).reference_charge * ATOM(j).rdists[k];
			}
			#undef ATOM
		}
}

/* Find atom type for a particular atom */
int get_atom_type_idx(const struct atom * const a) {

	assert(a != NULL);

	for(int i = 0; i < ts.atom_types_count; i++)
		if(at_compare_against_atom(&ts.atom_types[i], a))
			return i;

	return NOT_FOUND;
}

/* Do some preprocessing to simplify things later on */
void preprocess_molecules(void) {

	if(s.mode == MODE_PARAMS || s.mode == MODE_CROSS) {
		/* Calculate sum and average of the charges in the molecule */
		for(int i = 0; i < ts.molecules_count; i++)
			m_calculate_charge_stats(&ts.molecules[i]);
	}

	if(s.mode == MODE_PARAMS) {
		/* Calculate average electronegativies */
		for(int i = 0; i < ts.molecules_count; i++)
			m_calculate_avg_electronegativity(&ts.molecules[i]);

		/* Calculate reciprocal distances of atoms for all molecules */
		for(int i = 0; i < ts.molecules_count; i++)
			m_calculate_rdists(&ts.molecules[i]);

		/* Calculate auxiliary sum */
		calculate_y();
	}

	/* Finally, fill indices to atoms of a particular kind */
	fill_atom_types();
}

/* Prints information about the training set */
void ts_info(void) {

	printf("\nTraining set info\n");

	printf("Molecules: %5d  Atoms: %8d  Atom types: %2d\n", ts.molecules_count, ts.atoms_count, ts.atom_types_count);
	printf("Atom type     # atoms       %% atoms     # molecules\n");

	for(int i = 0; i < ts.atom_types_count; i++) {
		#define AT ts.atom_types[i]
		char buff[10];
		at_format_text(&AT, buff);
		printf(" %-10s %8d        %6.3f     %8d\n", buff, AT.atoms_count, 100.0f * (float) AT.atoms_count / ts.atoms_count, AT.molecules_count);
		#undef AT
	}

	printf("\n");
}

/* Format text caption for a particular atom type */
void at_format_text(const struct atom_type * const at, char * const buff) {

	assert(at != NULL);
	assert(buff != NULL);

	/* Note that buff should have size at least 9, this is not checked! */
	switch(s.at_customization) {
		case AT_CUSTOM_ELEMENT:
			snprintf(buff, 10, "%s", convert_Z_to_symbol(at->Z));
			break;
		case AT_CUSTOM_ELEMENT_BOND:
			snprintf(buff, 10, "%s %1d", convert_Z_to_symbol(at->Z), at->bond_order);
			break;
		case AT_CUSTOM_USER:
			snprintf(buff, 10, "%s", at->type_string);
			break;
		default:
			/* Something bad happened */
			assert(0);
	}
}

int get_atom_type_idx_from_text(const char * const str) {

	assert(str != NULL);

	struct atom a;

	char symbol[3];
	int bonds;

	switch(s.at_customization) {
		case AT_CUSTOM_ELEMENT:
			sscanf(str, "%2s\n", symbol);
			a.Z = convert_symbol_to_Z(symbol);
			break;
		case AT_CUSTOM_ELEMENT_BOND:
			sscanf(str, "%2s %d\n", symbol, &bonds);
			a.Z = convert_symbol_to_Z(symbol);
			a.bond_order = bonds;
			break;
		case AT_CUSTOM_USER:
			strncpy(a.type_string, str, 10);
			break;
		default:
			/* Something bad happened */
			assert(0);
	}

	/* Check if the conversion was succesful */
	if(a.Z == 0)
		EXIT_ERROR(RUN_ERROR, "Cannot convert \"%s\" to a atom type description!\n", str);

	return get_atom_type_idx(&a);
}

/* Fill necessary atom type items according to an atom */
static void at_fill_from_atom(struct atom_type * const at, const struct atom * const a) {

	assert(at != NULL);
	assert(a != NULL);

	switch(s.at_customization) {
			case AT_CUSTOM_ELEMENT:
				at->Z = a->Z;
				break;
			case AT_CUSTOM_ELEMENT_BOND:
				at->Z = a->Z;
				at->bond_order = a->bond_order;
				break;
			case AT_CUSTOM_USER:
				at->Z = a->Z;
				strncpy(at->type_string, a->type_string, 10);
				break;
			default:
				/* Something bad happened */
				assert(0);
	}
}

/* Compare if the provided atom is of a particular atom type */
static int at_compare_against_atom(const struct atom_type * const at, const struct atom * const a) {

	assert(at != NULL);
	assert(a != NULL);

	switch(s.at_customization) {
			case AT_CUSTOM_ELEMENT:
				return at->Z == a->Z;
			case AT_CUSTOM_ELEMENT_BOND:
				return at->Z ==  a->Z && at->bond_order == a->bond_order;
			case AT_CUSTOM_USER:
				return !strncmp(at->type_string, a->type_string, 10);
			default:
				/* Something bad happened */
				assert(0);
	}
}

/* List molecules for which we don't have parameters */
static void list_molecules_without_parameters(void) {

	/* Detect those molecules */
	for(int i = 0; i < ts.atom_types_count; i++)
		if(!ts.atom_types[i].has_parameters) {
			for(int j = 0; j < ts.atom_types[i].atoms_count; j++) {
				ts.molecules[ts.atom_types[i].atoms_molecule_idx[j]].has_parameters = 0;
			}
		}

	/* Print some statistics */
	int without_parameters_count = 0;
	for(int i = 0; i < ts.molecules_count; i++)
		if(!ts.molecules[i].has_parameters)
			without_parameters_count++;

	printf("\nLoaded parameters covering %d out of %d molecules (%4.2f %%).\n",
	ts.molecules_count - without_parameters_count, ts.molecules_count,
	100.0f * (ts.molecules_count - without_parameters_count) / ts.molecules_count);

	if(!without_parameters_count || !s.list_omitted_molecules)
		return;

	/* And the affected molecules themselves */
	printf("List of molecules without appropriate parameters loaded:\n");
	for(int i = 0; i < ts.molecules_count; i++)
		if(!ts.molecules[i].has_parameters)
			printf("%s; ", ts.molecules[i].name);

	printf("\n");
}

/* List molecules for which the charges haven't been loaded */
static void list_molecules_without_charges(void) {

	/* Print some statistics */
	int without_charges_count = 0;
	for(int i = 0; i < ts.molecules_count; i++)
		if(!ts.molecules[i].has_charges)
			without_charges_count++;

	printf("\nLoaded charges covering %d out of %d molecules (%4.2f %%).\n",
		ts.molecules_count - without_charges_count, ts.molecules_count,
		100.0f * (ts.molecules_count - without_charges_count) / ts.molecules_count);

	if(!without_charges_count || !s.list_omitted_molecules)
		return;

	/* And the affected molecules themselves */
	printf("List of molecules without appropriate charges loaded:\n");
	for(int i = 0; i < ts.molecules_count; i++)
		if(!ts.molecules[i].has_charges)
			printf("%s; ", ts.molecules[i].name);

	printf("\n");
}

/* List molecules which have wrong structure */
static void list_invalid_molecules(void) {

	/* Print some statistics */
	int invalid_molecules_count = 0;
	for(int i = 0; i < ts.molecules_count; i++)
		if(!ts.molecules[i].is_valid)
			invalid_molecules_count++;

	printf("\nLoaded %d valid molecules out of total %d (%4.2f %%).\n",
		ts.molecules_count - invalid_molecules_count, ts.molecules_count,
		100.0f * (ts.molecules_count - invalid_molecules_count) / ts.molecules_count);

	if(!invalid_molecules_count || !s.list_omitted_molecules)
		return;

	/* And the affected molecules themselves */
	printf("List of invalid molecules:\n");
	for(int i = 0; i < ts.molecules_count; i++)
		if(!ts.molecules[i].is_valid)
			printf("%s; ", ts.molecules[i].name);

	printf("\n");
}

/* List molecules which don't have right charges */
static void list_molecules_without_atom_types(void) {

	/* Print some statistics */
	int without_atom_types_count = 0;
	for(int i = 0; i < ts.molecules_count; i++)
		if(!ts.molecules[i].has_atom_types)
			without_atom_types_count++;

	printf("\nLoaded atom types for %d molecules out of total %d (%4.2f %%).\n",
		ts.molecules_count - without_atom_types_count, ts.molecules_count,
		100.0f * (ts.molecules_count - without_atom_types_count) / ts.molecules_count);

	if(!without_atom_types_count || !s.list_omitted_molecules)
		return;

	/* And the affected molecules themselves */
	printf("List of molecules without atom types:\n");
	for(int i = 0; i < ts.molecules_count; i++)
		if(!ts.molecules[i].has_atom_types)
			printf("%s; ", ts.molecules[i].name);

	printf("\n");
}

/* Remove molecules which can't be further used */
void discard_invalid_molecules_or_without_charges_or_parameters(void) {

	list_invalid_molecules();
	if(s.at_customization == AT_CUSTOM_USER)
		list_molecules_without_atom_types();

	if(s.mode == MODE_CHARGES || s.mode == MODE_CROSS || s.mode == MODE_COVER)
		list_molecules_without_parameters();

	if(s.mode == MODE_PARAMS || s.mode == MODE_CROSS)
		list_molecules_without_charges();

	/* Discard those molecules */
	int idx = 0;
	int number_of_discarded = 0;

	while(idx < ts.molecules_count) {
		int cond = 0;
		/* Different combination of parameters/charges are required for each mode, so act accordingly */
		if(s.mode == MODE_CHARGES)
			cond = !ts.molecules[idx].has_parameters || !ts.molecules[idx].is_valid;
		else if (s.mode == MODE_COVER)
			cond = !ts.molecules[idx].has_parameters;
		else if (s.mode == MODE_PARAMS)
			cond = !ts.molecules[idx].has_charges || !ts.molecules[idx].is_valid;
		else if (s.mode == MODE_CROSS)
			cond = !ts.molecules[idx].has_parameters || !ts.molecules[idx].has_charges || !ts.molecules[idx].is_valid;

		if(s.at_customization == AT_CUSTOM_USER)
			cond |= !ts.molecules[idx].has_atom_types;

		if(cond) {
			number_of_discarded++;
			/* Destroy molecule without parameters; fill its space with the last one */
			ts.atoms_count -= ts.molecules[idx].atoms_count;
			m_destroy(&ts.molecules[idx]);
			ts.molecules_count--;
			/* Check if the last molecule was discarded */
			if(idx != ts.molecules_count)
				ts.molecules[idx] = ts.molecules[ts.molecules_count];
		} else
			idx++;
	}

	if(number_of_discarded)
		printf("\nDiscarded %d molecules.\n", number_of_discarded);

	/* Shrink molecules array */
	ts.molecules = (struct molecule *) realloc(ts.molecules, sizeof(struct molecule) * ts.molecules_count);

	/* We need to rebuild atom types info */
	for(int i = 0; i < ts.atom_types_count; i++)
		at_destroy(&ts.atom_types[i]);

	free(ts.atom_types);

	fill_atom_types();
}

/* Print molecular formula to the string */
void get_sum_formula(const struct molecule * const m, char * buff, int n) {

	assert(m != NULL);
	assert(buff != NULL);

	int *atoms = (int *) calloc(103, sizeof(int));

	for(int i = 0; i < m->atoms_count; i++)
		atoms[m->atoms[i].Z] += 1;

	int len = 0;
	for(int i = 0; i < 103; i++) {
		if(atoms[i] == 1)
			len = snprintf(buff, n - len, "%s", convert_Z_to_symbol(i));
		else if(atoms[i] > 1)
			len = snprintf(buff, n - len, "%s%d", convert_Z_to_symbol(i), atoms[i]);
		else
			continue;

		buff += len;
	}

	free(atoms);
}
