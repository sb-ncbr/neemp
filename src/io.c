/*
 * NEEMP - io.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013, 2014
 *
 * */

#include <assert.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <math.h>
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

		if(i == MAX_MOLECULES)
			EXIT_ERROR(RUN_ERROR, "Maximum number of molecules (%d) reached. "
					      "Increase value of MAX_MOLECULES in config.h and recompile.\n", MAX_MOLECULES);
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
		if(idx == NOT_FOUND) {
			fprintf(stderr, "Molecule %s not loaded from .sdf file. Skipping the charge record (%s).\n", line, s.chg_file);
			/* Skip the whole record */
			do {
				if(!fgets(line, MAX_LINE_LEN, f)) {
					if(feof(f))
						break;
					else
						EXIT_ERROR(IO_ERROR, "Reading failed when skipping record (%s).\n", s.chg_file);
				}
			} while(strcmp(line, "\n"));

			/* Go to the next record */
			continue;
		}

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

static xmlNodePtr get_child_node_by_name(xmlNodePtr node, const char * const name) {

	xmlNodePtr curr_node = node->children;

	while(curr_node != NULL) {

		if(!strcmp((char *) curr_node->name, name))
			return curr_node;

		curr_node = curr_node->next;
	}

	return NULL;
}

void load_parameters(struct kappa_data * const kd) {

	assert(kd != NULL);

	xmlDocPtr doc = NULL;
	xmlNodePtr root_node = NULL;

	if((doc = xmlReadFile(s.par_file, NULL, XML_PARSE_NOBLANKS)) == NULL)
		EXIT_ERROR(IO_ERROR, "Cannot open or parse .par file \"%s\".\n", s.par_file);

	root_node = xmlDocGetRootElement(doc);


	xmlNodePtr parameters_node = get_child_node_by_name(root_node, "Parameters");
	if(parameters_node == NULL)
		EXIT_ERROR(IO_ERROR, "%s", "Ill-formed .par file. No Parameters node.\n");

	xmlChar *atom_type = xmlGetProp(parameters_node, BAD_CAST "AtomType");
	if(atom_type == NULL) {
		/* Assume ElemBond by default */
		atom_type = (xmlChar *) malloc(sizeof(xmlChar) * 10);
		snprintf((char *) atom_type, 9, "%s", "ElemBond");
	}

	/* Check if the command-line settings matches the entry in the .par file */
	if(strcmp((char *) atom_type, get_atom_types_by_string(s.at_customization)))
		EXIT_ERROR(RUN_ERROR, "atom-types-by \"%s\" doesn't match with provided settings \"%s\".\n",
			(char *) atom_type, get_atom_types_by_string(s.at_customization));

	xmlChar *kappa = xmlGetProp(parameters_node, BAD_CAST "Kappa");

	if(kappa == NULL)
		EXIT_ERROR(IO_ERROR, "%s", "Ill-formed .par file. No Kappa property.\n");

	kd->kappa = (float) atof((char *) kappa);

	xmlFree(kappa);
	xmlFree(atom_type);

	xmlNodePtr element_node = parameters_node->children;

	while(element_node != NULL) {

		xmlChar *symbol = xmlGetProp(element_node, BAD_CAST "Name");
		if(symbol == NULL)
			EXIT_ERROR(IO_ERROR, "%s", "Ill-formed .par file. No Name property.\n");

		xmlNodePtr ab_node = element_node->children;
		while(ab_node != NULL) {

			xmlChar *parameter_a = NULL, *parameter_b = NULL;

			char buff[10];
			switch(s.at_customization) {

				case AT_CUSTOM_ELEMENT:
					snprintf(buff, 10, "%2s", (char *) symbol);
					break;
				case AT_CUSTOM_ELEMENT_BOND: {
					xmlChar *bond_order = xmlGetProp(ab_node, BAD_CAST "Type");
					if(!bond_order)
						EXIT_ERROR(IO_ERROR, "Could not load parameters for element %s\n", (char *) symbol);

					int bond = atoi((char *) bond_order);

					snprintf(buff, 10, "%2s %1d", (char *) symbol, bond);
					xmlFree(bond_order);
					break;
				}
				default:
					assert(0);
			}

			parameter_a = xmlGetProp(ab_node, BAD_CAST "A");
			parameter_b = xmlGetProp(ab_node, BAD_CAST "B");

			if(!parameter_a || !parameter_b)
				EXIT_ERROR(IO_ERROR, "Could not load parameters for element %s\n", (char *) symbol);

			int atom_type_idx = get_atom_type_idx_from_text(buff);
			if(atom_type_idx != NOT_FOUND) {
				/* Store alpha and beta parameters */
				kd->parameters_alpha[atom_type_idx] = (float) atof((char *) parameter_a);
				kd->parameters_beta[atom_type_idx] = (float) atof((char *) parameter_b);
			}

			xmlFree(parameter_a);
			xmlFree(parameter_b);

			ab_node = ab_node->next;
		}

		xmlFree(symbol);

		element_node = element_node->next;
	}

	xmlFreeDoc(doc);
	xmlCleanupParser();

	/* Check if we load all necessary parameters */
	for(int i = 0; i < ts.atom_types_count; i++) {
		if(fabsf(kd->parameters_alpha[i] * kd->parameters_beta[i]) <= 0.0f) {
			char buff[10];
			at_format_text(&ts.atom_types[i], buff);
			EXIT_ERROR(RUN_ERROR, "No parameters loaded for: %s\n", buff);
		}
	}
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

	sscanf(line + 34, "%5s\n", version);

	if(!strcmp(version, "V2000")) {
		m->atoms_count = strn2int(line, 3);
		bonds_count = strn2int(line + 3, 3);

		/* Perform some checks on the values read */
		if(m->atoms_count > MAX_ATOMS_PER_MOLECULE)
			EXIT_ERROR(IO_ERROR, "Number of atoms is incorrect for molecule \"%s\" (%s).\n", m->name, s.sdf_file);

		if(bonds_count > MAX_BONDS_PER_MOLECULE)
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


			if(s.mode == MODE_PARAMS) {
				m->atoms[i].rdists = (double *) calloc(m->atoms_count, sizeof(double));
				if(!m->atoms[i].rdists)
					EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for atom distances.\n");
			}

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
		} while(strncmp(line, "$$$$", strlen("$$$$")));

	} else if(!strcmp(version, "V3000")) {

		/* Read BEGIN CTAB entry */
		if(!fgets(line, MAX_LINE_LEN, f))
			EXIT_ERROR(IO_ERROR, "Reading failed for the molecule \"%s\" (%s).\n", m->name, s.sdf_file);

		if(strncmp(line, "M  V30 BEGIN CTAB", strlen("M  V30 BEGIN CTAB")))
			EXIT_ERROR(IO_ERROR, "Incorrect format of the BEGIN CTAB entry for the molecule \"%s\" (%s).\n", m->name, s.sdf_file);

		/* Read COUNTS entry */
		if(!fgets(line, MAX_LINE_LEN, f))
			EXIT_ERROR(IO_ERROR, "Reading failed for the molecule \"%s\" (%s).\n", m->name, s.sdf_file);

		if(strncmp(line, "M  V30 COUNTS", strlen("M  V30 COUNTS")))
			EXIT_ERROR(IO_ERROR, "Incorrect format of the COUNTS entry for the molecule \"%s\" (%s).\n", m->name, s.sdf_file);

		sscanf(line + 14, "%d %d", &m->atoms_count, &bonds_count);

		/* Perform some checks on the values read */
		if(m->atoms_count > MAX_ATOMS_PER_MOLECULE)
			EXIT_ERROR(IO_ERROR, "Number of atoms is incorrect for molecule \"%s\" (%s).\n", m->name, s.sdf_file);

		if(bonds_count > MAX_BONDS_PER_MOLECULE)
			EXIT_ERROR(IO_ERROR, "Number of bonds is incorrect for molecule \"%s\" (%s).\n", m->name, s.sdf_file);

		m->atoms = (struct atom *) malloc(sizeof(struct atom) * m->atoms_count);
		if(!m->atoms)
			EXIT_ERROR(MEM_ERROR, "Cannot allocate memory for atoms in molecule \"%s\" (%s).\n", m->name, s.sdf_file);

		/* Read BEGIN ATOM entry */
		if(!fgets(line, MAX_LINE_LEN, f))
			EXIT_ERROR(IO_ERROR, "Reading failed for the molecule \"%s\" (%s).\n", m->name, s.sdf_file);

		if(strncmp(line, "M  V30 BEGIN ATOM", strlen("M  V30 BEGIN ATOM")))
			EXIT_ERROR(IO_ERROR, "Incorrect format of the BEGIN ATOM entry for the molecule \"%s\" (%s).\n", m->name, s.sdf_file);

		/* Process individual atom records */
		for(int i = 0; i < m->atoms_count; i++) {
			if(!fgets(line, MAX_LINE_LEN, f))
				EXIT_ERROR(IO_ERROR, "Reading failed for the molecule \"%s\" (%s).\n", m->name, s.sdf_file);

			char atom_symbol[3];
			int tmp;

			if(strncmp(line, "M  V30", strlen("M  V30")))
				EXIT_ERROR(IO_ERROR, "Incorrect format of ATOM entry for the molecule \"%s\" (%s).\n", m->name, s.sdf_file);

			sscanf(line + 7, "%d %2s %f %f %f", &tmp, atom_symbol,\
				&m->atoms[i].position[0], &m->atoms[i].position[1], &m->atoms[i].position[2]);

			if(s.mode == MODE_PARAMS) {
				m->atoms[i].rdists = (double *) calloc(m->atoms_count, sizeof(double));
				if(!m->atoms[i].rdists)
					EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for atom distances.\n");
			}

			m->atoms[i].Z = convert_symbol_to_Z(atom_symbol);
			if(m->atoms[i].Z == 0)
				EXIT_ERROR(IO_ERROR, "Invalid element \"%s\" in the molecule \"%s\" (%s).\n", atom_symbol, m->name, s.sdf_file);

			m->atoms[i].bond_order = 1;
		}

		/* Read END ATOM entry */
		if(!fgets(line, MAX_LINE_LEN, f))
			EXIT_ERROR(IO_ERROR, "Reading failed for the molecule \"%s\" (%s).\n", m->name, s.sdf_file);

		if(strncmp(line, "M  V30 END ATOM", strlen("M  V30 END ATOM")))
			EXIT_ERROR(IO_ERROR, "Incorrect format of the END ATOM entry for the molecule \"%s\" (%s).\n", m->name, s.sdf_file);

		/* Read BEGIN BOND entry */
		if(!fgets(line, MAX_LINE_LEN, f))
			EXIT_ERROR(IO_ERROR, "Reading failed for the molecule \"%s\" (%s).\n", m->name, s.sdf_file);

		if(strncmp(line, "M  V30 BEGIN BOND", strlen("M  V30 BEGIN BOND")))
			EXIT_ERROR(IO_ERROR, "Incorrect format of the BEGIN BOND entry for the molecule \"%s\" (%s).\n", m->name, s.sdf_file);

		for(int i = 0; i < bonds_count; i++) {
			if(!fgets(line, MAX_LINE_LEN, f))
				EXIT_ERROR(IO_ERROR, "Reading failed for the molecule \"%s\" (%s).\n", m->name, s.sdf_file);

			int tmp, bond_order, atom1, atom2;

			if(strncmp(line, "M  V30", strlen("M  V30")))
				EXIT_ERROR(IO_ERROR, "Incorrect format of ATOM entry for the molecule \"%s\" (%s).\n", m->name, s.sdf_file);

			sscanf(line + 7, "%d %d %d %d", &tmp, &bond_order, &atom1, &atom2);

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

		/* Read END BOND entry */
		if(!fgets(line, MAX_LINE_LEN, f))
			EXIT_ERROR(IO_ERROR, "Reading failed for the molecule \"%s\" (%s).\n", m->name, s.sdf_file);

		if(strncmp(line, "M  V30 END BOND", strlen("M  V30 END BOND")))
			EXIT_ERROR(IO_ERROR, "Incorrect format of the END BOND entry for the molecule \"%s\" (%s).\n", m->name, s.sdf_file);

		/* Skip other entries */
		do {
			if(!fgets(line, MAX_LINE_LEN, f))
				EXIT_ERROR(IO_ERROR, "Reading failed for the molecule \"%s\" (%s).\n", m->name, s.sdf_file);
		} while(strncmp(line, "M  V30 END CTAB", strlen("M  V30 END CTAB")));

		/* Skip the rest of the record */
		do {
			if(!fgets(line, MAX_LINE_LEN, f))
				EXIT_ERROR(IO_ERROR, "Reading failed for the molecule \"%s\" (%s).\n", m->name, s.sdf_file);
		} while(strncmp(line, "$$$$", strlen("$$$$")));
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
		fprintf(f, "R: %6.4f  RMSD: %4.2e  MSE: %4.2e  D_avg: %4.2e  D_max: %4.2e\n",
			ss->best->per_molecule_stats[i].R, ss->best->per_molecule_stats[i].RMSD, ss->best->per_molecule_stats[i].MSE,
			ss->best->per_molecule_stats[i].D_avg, ss->best->per_molecule_stats[i].D_max);
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

	xmlDocPtr doc = NULL;
	xmlNodePtr root_node = NULL;

	doc = xmlNewDoc(BAD_CAST "1.0");

	root_node = xmlNewNode(NULL, BAD_CAST "ParameterSet");
	xmlDocSetRootElement(doc, root_node);

	xmlNodePtr params_node = xmlNewChild(root_node, NULL, BAD_CAST "Parameters", NULL);

	char buff[10];
	snprintf(buff, 10, "%s", get_atom_types_by_string(s.at_customization));

	xmlNewProp(params_node, BAD_CAST "AtomType", BAD_CAST buff);

	snprintf(buff, 10, "%6.4f", ss->best->kappa);

	xmlNewProp(params_node, BAD_CAST "Kappa", BAD_CAST buff);
	xmlAddChild(root_node, params_node);

	for(int i = 0; i < ts.atom_types_count; i++) {

		snprintf(buff, 10, "%2s", convert_Z_to_symbol(ts.atom_types[i].Z));

		xmlNodePtr element_node = NULL;
		xmlNodePtr curr_node = params_node->children;

		while(curr_node != NULL) {
			if(!strncmp((char *) xmlGetProp(curr_node, BAD_CAST "Name"), buff, 10))
				element_node = curr_node;

			curr_node = curr_node->next;
		}

		if(!element_node) {
			element_node = xmlNewChild(params_node, NULL, BAD_CAST "Element", NULL);
			xmlNewProp(element_node, BAD_CAST "Name", BAD_CAST buff);
		}

		xmlNodePtr bond_node = xmlNewChild(element_node, NULL, BAD_CAST "Bond", NULL);

		if(s.at_customization == AT_CUSTOM_ELEMENT_BOND) {
			snprintf(buff, 10, "%d", ts.atom_types[i].bond_order);
			xmlNewProp(bond_node, BAD_CAST "Type", BAD_CAST buff);
		}

		snprintf(buff, 10, "%6.4f", ss->best->parameters_alpha[i]);
		xmlNewProp(bond_node, BAD_CAST "A", BAD_CAST buff);
		snprintf(buff, 10,"%6.4f", ss->best->parameters_beta[i]);
		xmlNewProp(bond_node, BAD_CAST "B", BAD_CAST buff);
	}

	if(xmlSaveFormatFile(s.par_out_file, doc, 1) == -1)
		EXIT_ERROR(IO_ERROR, "Cannot open file %s for writing the parameters.\n", s.par_out_file);

	xmlFreeDoc(doc);
}
