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

	FILE * const f = fopen(s.sdf_filename, "r");
	if(!f)
		EXIT_ERROR(IO_ERROR, "Cannot open .sdf file \"%s\".\n", s.sdf_filename);

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

	FILE * const f = fopen(s.chg_filename, "r");
	if(!f)
		EXIT_ERROR(IO_ERROR, "Cannot open .chg file \"%s\".\n", s.chg_filename);

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
		fgets(line, MAX_LINE_LEN, f);
		sscanf(line, "%d", &atoms_count);
		if(atoms_count != ts.molecules[idx].atoms_count)
			EXIT_ERROR(IO_ERROR, "Number of atoms in molecule \"%s\" don't match.\n", ts.molecules[idx].name);

		/* Load actual charges */
		for(int i = 0; i < atoms_count; i++) {
			fgets(line, MAX_LINE_LEN, f);
			int tmp_int;
			char tmp_str[2];
			sscanf(line, "%d %s %f\n", &tmp_int, tmp_str, &ts.molecules[idx].atoms[i].reference_charge);
		}

		ts.molecules[idx].charges_loaded = 1;

		/* Read empty line */
		fgets(line, MAX_LINE_LEN, f);
	}
	fclose(f);
}

/* TODO */
void load_parameters(void) {

	FILE * const f = fopen(s.par_filename, "r");
	if(!f)
		EXIT_ERROR(IO_ERROR, "Cannot open .par file \"%s\".\n", s.par_filename);

	fclose(f);
}

/* Convert n characters of a string to int */
static int strn2int(const char * const str, int n) {

	char buff[n];
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
	fgets(line, MAX_LINE_LEN, f);
	/* 3rd line is for comments, skip it */
	fgets(line, MAX_LINE_LEN, f);

	/* Read Counts Line
	 *
	 * format: aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
	 * aaa - number of atoms
	 * bbb - number of bonds
	 * vvvvvv - version (either V2000 or V3000)
	 * the rest is not used by NEEMP */

	int bonds_count;
	char version[MAX_LINE_LEN];

	fgets(line, MAX_LINE_LEN, f);
	m->atoms_count = strn2int(line, 3);
	bonds_count = strn2int(line + 3, 3);
	sscanf(line + 33, "%6s\n", version);

	if(!strcmp(version, "V2000")) {
		/* Perform some checks on the values read */
		if(MIN_ATOMS_PER_MOLECULE > m->atoms_count || m->atoms_count > MAX_ATOMS_PER_MOLECULE)
			EXIT_ERROR(IO_ERROR, "Number of atoms is incorrect for molecule \"%s\"\n", m->name);

		if(MIN_BONDS_PER_MOLECULE > bonds_count || bonds_count > MAX_BONDS_PER_MOLECULE)
			EXIT_ERROR(IO_ERROR, "Number of bonds is incorrect for molecule \"%s\"\n", m->name);

		m->atoms = (struct atom *) malloc(sizeof(struct atom) * m->atoms_count);
		if(!m->atoms)
			EXIT_ERROR(MEM_ERROR, "Cannot allocate memory for atoms in molecule \"%s\".\n", m->name);

		/* Process Atom Block
		 *
		 * format: xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
		 * x, y, z - coordinates
		 * aaa - atom symbol
		 * the rest is not used by NEEMP */

		for(int i = 0; i < m->atoms_count; i++) {
			char atom_symbol[3];
			fgets(line, MAX_LINE_LEN, f);
			sscanf(line, "%f %f %f %s", &m->atoms[i].position[0], &m->atoms[i].position[1], &m->atoms[i].position[2], atom_symbol);

			m->atoms[i].rdists = (double *) calloc(m->atoms_count, sizeof(double));
			if(!m->atoms[i].rdists)
				EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for atom distances.\n");

			m->atoms[i].Z = convert_symbol_to_Z(atom_symbol);
			if(m->atoms[i].Z == 0)
				EXIT_ERROR(IO_ERROR, "Invalid element \"%s\" in the molecule \"%s\".\n", atom_symbol, m->name);

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

			fgets(line, MAX_LINE_LEN, f);

			atom1 = strn2int(line, 3);
			atom2 = strn2int(line + 3, 3);
			bond_order = strn2int(line + 6, 3);

			/* Perform some checks on the data */
			if(atom1 > m->atoms_count || atom2 > m->atoms_count)
				EXIT_ERROR(IO_ERROR, "Invalid atom number in the molecule \"%s\".\n", m->name);

			if(bond_order > 3)
				EXIT_ERROR(IO_ERROR, "Invalid bond order in the molecule \"%s\".\n", m->name);

			/* Adjust bond orders of the atoms */
			if(m->atoms[atom1 - 1].bond_order < bond_order)
				m->atoms[atom1 -1].bond_order = bond_order;

			if(m->atoms[atom2 - 1].bond_order < bond_order)
				m->atoms[atom2 -1].bond_order = bond_order;
		}

		/* Skip rest of the record */
		do {
			fgets(line, MAX_LINE_LEN, f);
		} while(strncmp(line, "$$$$", 4));

	} else if(!strcmp(version, "V3000")) {
		/* TODO */
		EXIT_ERROR(IO_ERROR, "%s", "MDL files with V3000 format are currently unsupported.\n");
	} else
		EXIT_ERROR(IO_ERROR, "MDL file with unknown version \"%s\".\n", version);

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


/* Output reference charges, EEM charges and their differences */
void output_charges_stats(const struct subset * const ss) {

	assert(ss != NULL);
	assert(ss->best != NULL);

	FILE *f = fopen(s.chgout_filename, "w");
	if(!f)
		EXIT_ERROR(IO_ERROR, "Cannot open file %s for writing the charges stats.\n", s.chgout_filename);

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
