/*
 * NEEMP - neemp.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#include <string.h>

#include "neemp.h"
#include "parser.h"
#include "settings.h"
#include "structures.h"

struct training_set ts;
struct settings s;

int main(int argc, char **argv) {

	strcpy(s.sdf_filename, "../molecules.sdf");
	strcpy(s.chg_filename, "../B3LYP_STO-3G.mchrg");

	load_molecules();
	load_charges();

	preprocess_molecules();

	ts_destroy();

	return RETURN_OK;
}
