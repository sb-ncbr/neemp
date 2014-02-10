/*
 * NEEMP - io.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "neemp.h"
#include "io.h"
#include "settings.h"
#include "subset.h"
#include "structures.h"

extern const struct settings s;
extern struct training_set ts;

static int load_molecule(FILE * const f, struct molecule * const m);
static int find_molecule_by_name(const char * const name);
static int strn2int(const char * const str, int n);

/* Load all molecules from .sdf file */
void load_molecules(void) {

	FILE * const f = fopen(s.sdf_file, "r");
	if(!f)
		EXIT_ERROR(IO_ERROR, "Cannot open .sdf file \"%s\".\n", s.sdf_file);

	ts.molecules = (struct molecule *) malloc(sizeof(struct molecule) * MAX_MOLECULES);
	if(!ts.molecules)
		EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for molecules\n");

	/* Load molecules one by one */
	int i = 0;
	while(!load_molecule(f, &ts.molecules[i])) {
		ts.atoms_count += ts.molecules[i].atoms_count;
		i++;
	}

	fclose(f);

	/* Free unused memory */
	ts.molecules_count = i;
	ts.molecules = (struct molecule *) realloc(ts.molecules, sizeof(struct molecule) * ts.molecules_count);
}

/* Load atomic charges from .chg file */
void load_charges(void) {

	/* .chg file
	 *
	 * format:
	 * -NAME-OF-THE-MOLECULE-
	 * -NUMBER-OF-ATOMS-
	 * -ATOM-NUMBER- -ATOM-SYMBOL- -CHARGE-
	 * -ATOM-NUMBER- -ATOM-SYMBOL- -CHARGE-
	 * -ATOM-NUMBER- -ATOM-SYMBOL- -CHARGE-
	 * [etc.]
	 *
	 * [EMPTY LINE terminates the record]
	 */

	FILE * const f = fopen(s.chg_file, "r");
	if(!f)
		EXIT_ERROR(IO_ERROR, "Cannot open .chg file \"%s\".\n", s.chg_file);

	char line[MAX_LINE_LEN];
	memset(line, 0x0, MAX_LINE_LEN * sizeof(char));

	while(1) {
		/* Break if no data is available */
		if(!fgets(line, MAX_LINE_LEN, f))
			break;

		/* First line is the name; strip newline character */
		int len = strlen(line);
		line[len - 1] = '\0';

		/* Find corresponding previously loaded molecule */
		int idx = find_molecule_by_name(line);
		if(idx == NOT_FOUND)
			EXIT_ERROR(IO_ERROR, "Molecule \"%s\" was not loaded.\n", line);

		/* Check if numbers of atoms match*/
		int atoms_count;
		if(!fgets(line, MAX_LINE_LEN, f))
			EXIT_ERROR(IO_ERROR, "Reading failed for the molecule \"%s\" (%s).\n", ts.molecules[idx].name, s.chg_file);

		sscanf(line, "%d", &atoms_count);
		if(atoms_count != ts.molecules[idx].atoms_count)
			EXIT_ERROR(IO_ERROR, "Number of atoms in molecule \"%s\" don't match (%s).\n", ts.molecules[idx].name, s.chg_file);

		/* Load actual charges */
		for(int i = 0; i < atoms_count; i++) {
			if(!fgets(line, MAX_LINE_LEN, f))
				EXIT_ERROR(IO_ERROR, "Reading charges failed for the molecule \"%s\" (%s).\n", ts.molecules[idx].name, s.chg_file);

			int tmp_int;
			char tmp_str[2];
			sscanf(line, "%d %s %f\n", &tmp_int, tmp_str, &ts.molecules[idx].atoms[i].reference_charge);
		}

		ts.molecules[idx].charges_loaded = 1;

		/* Read empty line */
		if(!fgets(line, MAX_LINE_LEN, f) && !feof(f))
			EXIT_ERROR(IO_ERROR, "Reading empty separator line failed after the molecule \"%s\" (%s).\n", ts.molecules[idx].name, s.chg_file);
	}
	fclose(f);
}

void load_parameters(struct kappa_data * const kd) {

	assert(kd != NULL);

	FILE * const f = fopen(s.par_file, "r");
	if(!f)
		EXIT_ERROR(IO_ERROR, "Cannot open .par file \"%s\".\n", s.par_file);

	char line[MAX_LINE_LEN];
	memset(line, 0x0, MAX_LINE_LEN);

	/* Skip comments */
	do {
		if(!fgets(line, MAX_LINE_LEN, f))
			EXIT_ERROR(IO_ERROR, "Invalid format of \"%s\".\n", s.par_file);
	} while(line[0] == '#');

	sscanf(line, "%f\n", &kd->kappa);

	/* Read actual parameters */
	for(int i = 0; i < ts.atom_types_count; i++) {

		if(!fgets(line, MAX_LINE_LEN, f))
			EXIT_ERROR(IO_ERROR, "Invalid format of \"%s\".\n", s.par_file);

		/* Determine atom type from string */
		int atom_type_idx = get_atom_type_idx_from_text(line);

		if(atom_type_idx == NOT_FOUND)
			EXIT_ERROR(IO_ERROR, "Invalid format of atom type specifier on line \"%s\" (%s).", line, s.par_file);

		/* Read alpha and beta */
		sscanf(line + 9, "%f %f\n", &kd->parameters_alpha[atom_type_idx], &kd->parameters_beta[atom_type_idx]);
	}

	fclose(f);
}

/* Convert n characters of a string to int */
static int strn2int(const char * const str, int n) {

	char buff[n + 1];
	memset(buff, 0x0, n + 1);
	for(int i = 0; i < n; i++)
		buff[i] = str[i];

	return atoi(buff);
}

/* Load one molecule from a .sd file */
static int load_molecule(FILE * const f, struct molecule * const m) {

	assert(f != NULL);
	assert(m != NULL);

	/* Each molecule is stored in MOL format;
	 * for reference, see http://c4.cabrillo.edu/404/ctfile.pdf */

	char line[MAX_LINE_LEN];
	memset(line, 0x0, MAX_LINE_LEN * sizeof(char));

	/* Process 3-line Header Block */

	/* Do we reached EOF? */
	if(!fgets(line, MAX_LINE_LEN, f))
		return 1;

	/* 1st line is the name of the molecule */
	int len = strlen(line);
	m->name = (char *) calloc(len, sizeof(char));
	if(!m->name)
		EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for molecule name.\n");

	strncpy(m->name, line, len - 1);
	m->name[len - 1] = '\0';

	/* 2nd line contains some additional information, skip it */
	if(!fgets(line, MAX_LINE_LEN, f))
			EXIT_ERROR(IO_ERROR, "Reading failed for 2nd line of the molecule \"%s\" (%s).\n", m->name, s.sdf_file);
	/* 3rd line is for comments, skip it */
	if(!fgets(line, MAX_LINE_LEN, f))
			EXIT_ERROR(IO_ERROR, "Reading failed for 3rd line of the molecule \"%s\" (%s).\n", m->name, s.sdf_file);

	/* Read Counts Line
	 *
	 * format: aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
	 * aaa - number of atoms
	 * bbb - number of bonds
	 * vvvvvv - version (either V2000 or V3000)
	 * the rest is not used by NEEMP */

	int bonds_count;
	char version[MAX_LINE_LEN];

	if(!fgets(line, MAX_LINE_LEN, f))
			EXIT_ERROR(IO_ERROR, "Reading failed for 4th line of the molecule \"%s\" (%s).\n", m->name, s.sdf_file);
	m->atoms_count = strn2int(line, 3);
	bonds_count = strn2int(line + 3, 3);
	sscanf(line + 33, "%6s\n", version);

	if(!strcmp(version, "V2000")) {
		/* Perform some checks on the values read */
		if(MIN_ATOMS_PER_MOLECULE > m->atoms_count || m->atoms_count > MAX_ATOMS_PER_MOLECULE)
			EXIT_ERROR(IO_ERROR, "Number of atoms is incorrect for molecule \"%s\" (%s).\n", m->name, s.sdf_file);

		if(MIN_BONDS_PER_MOLECULE > bonds_count || bonds_count > MAX_BONDS_PER_MOLECULE)
			EXIT_ERROR(IO_ERROR, "Number of bonds is incorrect for molecule \"%s\" (%s).\n", m->name, s.sdf_file);

		m->atoms = (struct atom *) malloc(sizeof(struct atom) * m->atoms_count);
		if(!m->atoms)
			EXIT_ERROR(MEM_ERROR, "Cannot allocate memory for atoms in molecule \"%s\" (%s).\n", m->name, s.sdf_file);

		/* Process Atom Block
		 *
		 * format: xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
		 * x, y, z - coordinates
		 * aaa - atom symbol
		 * the rest is not used by NEEMP */

		for(int i = 0; i < m->atoms_count; i++) {
			char atom_symbol[3];
			if(!fgets(line, MAX_LINE_LEN, f))
					EXIT_ERROR(IO_ERROR, "Reading failed for atom %d in the molecule \"%s\" (%s).\n", i + 1, m->name, s.sdf_file);
			sscanf(line, "%f %f %f %s", &m->atoms[i].position[0], &m->atoms[i].position[1], &m->atoms[i].position[2], atom_symbol);

			m->atoms[i].rdists = (double *) calloc(m->atoms_count, sizeof(double));
			if(!m->atoms[i].rdists)
				EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for atom distances.\n");

			m->atoms[i].Z = convert_symbol_to_Z(atom_symbol);
			if(m->atoms[i].Z == 0)
				EXIT_ERROR(IO_ERROR, "Invalid element \"%s\" in the molecule \"%s\" (%s).\n", atom_symbol, m->name, s.sdf_file);

			m->atoms[i].bond_order = 1;
		}

		/* Process Bond Block
		 *
		 * format: 111222tttsssxxxrrrccc
		 * 111 - first atom number
		 * 222 - second atom number
		 * ttt - bond type (either 1, 2 or 3)
		 * the rest is not used by NEEMP */

		for(int i = 0; i < bonds_count; i++) {
			int atom1, atom2, bond_order;

			if(!fgets(line, MAX_LINE_LEN, f))
					EXIT_ERROR(IO_ERROR, "Reading failed for bond no. %d in the molecule \"%s\" (%s).\n", i + 1, m->name, s.sdf_file);

			atom1 = strn2int(line, 3);
			atom2 = strn2int(line + 3, 3);
			bond_order = strn2int(line + 6, 3);

			/* Perform some checks on the data */
			if(atom1 > m->atoms_count || atom2 > m->atoms_count)
				EXIT_ERROR(IO_ERROR, "Invalid atom number in the molecule \"%s\" (%s).\n", m->name, s.sdf_file);

			if(bond_order > 3)
				EXIT_ERROR(IO_ERROR, "Invalid bond order in the molecule \"%s\" (%s).\n", m->name, s.sdf_file);

			/* Adjust bond orders of the atoms */
			if(m->atoms[atom1 - 1].bond_order < bond_order)
				m->atoms[atom1 -1].bond_order = bond_order;

			if(m->atoms[atom2 - 1].bond_order < bond_order)
				m->atoms[atom2 -1].bond_order = bond_order;
		}

		/* Skip rest of the record */
		do {
			if(!fgets(line, MAX_LINE_LEN, f))
					EXIT_ERROR(IO_ERROR, "Reading failed for the molecule \"%s\" (%s).\n", m->name, s.sdf_file);
		} while(strncmp(line, "$$$$", 4));

	} else if(!strcmp(version, "V3000")) {
		/* TODO */
		EXIT_ERROR(IO_ERROR, "MDL files with V3000 format are currently unsupported (%s).\n", s.sdf_file);
	} else
		EXIT_ERROR(IO_ERROR, "MDL file with unknown version \"%s\" (%s).\n", version, s.sdf_file);

	m->charges_loaded = 0;

	return 0;
}

/* Find index of the molecule with given name */
static int find_molecule_by_name(const char * const name) {

	assert(name != NULL);

	for(int i = 0; i < ts.molecules_count; i++)
		if(!strcmp(name, ts.molecules[i].name))
			return i;

	/* Not found */
	return NOT_FOUND;
}

void output_charges(const struct subset * const ss) {

	assert(ss != NULL);
	assert(ss->best != NULL);

	FILE *f = fopen(s.chg_out_file, "w");
	if(!f)
		EXIT_ERROR(IO_ERROR, "Cannot open file %s for writing the charges stats.\n", s.chg_out_file);

	int atoms_processed = 0;
	for(int i = 0; i < ts.molecules_count; i++) {
		fprintf(f, "%s\n", ts.molecules[i].name);
		fprintf(f, "%d\n", ts.molecules[i].atoms_count);

		for(int j = 0; j < ts.molecules[i].atoms_count; j++) {
			#define ATOM ts.molecules[i].atoms[j]
			fprintf(f, "%4d\t%2s\t%9.6f\n", j + 1, convert_Z_to_symbol(ATOM.Z),\
				ss->best->charges[atoms_processed + j]);
			#undef ATOM
		}
		atoms_processed += ts.molecules[i].atoms_count;
		fprintf(f, "\n");
	}

	fclose(f);
}

/* Output reference charges, EEM charges and their differences */
void output_charges_stats(const struct subset * const ss) {

	assert(ss != NULL);
	assert(ss->best != NULL);

	FILE *f = fopen(s.chg_stats_out_file, "w");
	if(!f)
		EXIT_ERROR(IO_ERROR, "Cannot open file %s for writing the charges stats.\n", s.chg_stats_out_file);

	fprintf(f, "IDX      TYPE        A.I.             EEM            DIFF\n");

	int atoms_processed = 0;
	for(int i = 0; i < ts.molecules_count; i++) {
		fprintf(f, "\n");
		fprintf(f, "%s\n", ts.molecules[i].name);
		fprintf(f, "%d\n", ts.molecules[i].atoms_count);
		for(int j = 0; j < ts.molecules[i].atoms_count; j++) {
			#define ATOM ts.molecules[i].atoms[j]
			fprintf(f, "%4d\t%2s %1d\t%9.6f\t%9.6f\t%9.6f\n", j + 1, convert_Z_to_symbol(ATOM.Z), ATOM.bond_order,
				ATOM.reference_charge, ss->best->charges[atoms_processed + j], ATOM.reference_charge - ss->best->charges[atoms_processed + j]);
			#undef ATOM
		}
		atoms_processed += ts.molecules[i].atoms_count;
	}

	fclose(f);
}

void output_parameters(const struct subset * const ss) {

	assert(ss != NULL);
	assert(ss->best != NULL);

	FILE *f = fopen(s.par_out_file, "w");
	if(!f)
		EXIT_ERROR(IO_ERROR, "Cannot open file %s for writing the parameters.\n", s.par_out_file);


	fprintf(f, "# NEEMP parameters file\n");
	fprintf(f, "# Customization (--at-customization parameter value): ");
	switch(s.at_customization) {
		case AT_CUSTOM_ELEMENT:
			fprintf(f, "element\n");
			break;
		case AT_CUSTOM_ELEMENT_BOND:
			fprintf(f, "element_bond\n");
			break;
		case AT_CUSTOM_PARTNER:
			fprintf(f, "partner\n");
			break;
		case AT_CUSTOM_VALENCE:
			fprintf(f, "valence\n");
			break;
		default:
			/* We should not be here */
			assert(0);
	}

	fprintf(f, "%6.4f\n", ss->best->kappa);

	for(int i = 0; i < ts.atom_types_count; i++) {
		char buff[10];
		at_format_text(&ts.atom_types[i], buff);
		fprintf(f, "%s\t%7.4f\t%7.4f\n", buff, ss->best->parameters_alpha[i], ss->best->parameters_beta[i]);
	}

	fclose(f);
}
