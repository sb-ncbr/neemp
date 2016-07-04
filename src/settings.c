/* Copyright 2013-2016 Tomas Racek (tom@krab1k.net)
 *
 * This file is part of NEEMP.
 *
 * NEEMP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * NEEMP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with NEEMP. If not, see <http://www.gnu.org/licenses/>.
 */

#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <string.h>

#include "limits.h"
#include "neemp.h"
#include "settings.h"
#include "../externals/lhs/latin_random.h"

extern struct settings s;

static void print_help(void);
static void print_version(void);

static char *atom_types_by_strings[] = {"Element", "ElemBond", "User"};

static struct option long_options[] = {

	{"help", no_argument, 0, 'h'},
	{"mode", required_argument, 0, 'm'},
	{"params-method", required_argument, 0, 'p'},
	{"verbose", no_argument, 0, 'v'},
	{"discard", required_argument, 0, 'd'},
	{"sort-by", required_argument, 0, 's'},
	{"version", no_argument, 0, 129},
	{"sdf-file", required_argument, 0, 130},
	{"par-file", required_argument, 0, 131},
	{"chg-file", required_argument, 0, 132},
	{"chg-out-file", required_argument, 0, 133},
	{"chg-stats-out-file", required_argument, 0, 134},
	{"par-out-file", required_argument, 0, 135},
	{"atb-file", required_argument, 0, 136},
	{"random-seed", required_argument, 0, 137},
	{"kappa-max", required_argument, 0, 140},
	{"kappa", required_argument, 0, 141},
	{"kappa-preset", required_argument, 0, 142},
	{"fs-precision", required_argument, 0, 143},
	{"atom-types-by", required_argument, 0, 151},
	{"tabu-size", required_argument, 0, 152},
	{"limit-iters", required_argument, 0, 160},
	{"limit-time", required_argument, 0, 161},
	{"check-charges", no_argument, 0, 170},
	{"max-threads", required_argument, 0, 171},
	{"list-omitted-molecules", no_argument, 0, 172},
	{"extra-precise", no_argument, 0, 173},
	{"om-pop-size", required_argument, 0, 180},
	{"de-f", required_argument, 0, 181},
	{"de-cr", required_argument, 0, 182},
	{"om-iters-max", required_argument, 0, 183},
	{"de-dither", no_argument, 0, 186},
	{"om-fix-kappa",required_argument, 0, 188},
	{"om-threads",required_argument, 0, 189},
	{"om-polish", required_argument, 0, 190},
	{"gm-iterations-beg", required_argument, 0, 192},
	{"gm-iterations-end", required_argument, 0, 193},
	{NULL, 0, 0, 0}
};


char *get_atom_types_by_string(enum atom_type_customization atc) {

	return atom_types_by_strings[atc];
}

/* Initialize default settings */
void s_init(void) {

	memset(s.sdf_file, 0x0, MAX_PATH_LEN * sizeof(char));
	memset(s.chg_file, 0x0, MAX_PATH_LEN * sizeof(char));
	memset(s.par_file, 0x0, MAX_PATH_LEN * sizeof(char));
	memset(s.atb_file, 0x0, MAX_PATH_LEN * sizeof(char));
	memset(s.par_out_file, 0x0, MAX_PATH_LEN * sizeof(char));
	memset(s.chg_out_file, 0x0, MAX_PATH_LEN * sizeof(char));
	memset(s.chg_stats_out_file, 0x0, MAX_PATH_LEN * sizeof(char));

	s.random_seed = -1;
	s.mode = MODE_NOT_SET;
	s.params_method = PARAMS_NOT_SET;
	s.full_scan_precision = 0.0f;
	s.kappa_max = 1.0f;
	s.full_scan_precision = 0.05f;
	s.population_size = 0;
	s.recombination_constant = -1;
	s.mutation_constant = -1;
	s.dither = 0;
	s.fixed_kappa = -1;
	s.om_threads = 1;
	s.om_iters = NO_LIMIT_ITERS;
	s.om_time = NO_LIMIT_TIME;
	s.polish = -1; /* 0 off, 1 only result, 2 result + during evolve, 3 result, evolve and some structures in initial population */
	s.gm_iterations_beg = 500;
	s.gm_iterations_end = 500;
	s.sort_by = SORT_NOT_SET;
	s.at_customization = AT_CUSTOM_ELEMENT_BOND;
	s.discard = DISCARD_OFF;
	s.tabu_size = 0.0f;
	s.limit_iters = NO_LIMIT_ITERS;
	s.limit_time = NO_LIMIT_TIME;
	s.check_charges = 0;
	s.max_threads = 1;
	s.list_omitted_molecules = 0;
	s.extra_precise = 0;
}

/* Prints help if --version is issued */
static void print_version(void) {

	printf("%s %s ", APP_NAME, APP_VERSION);
	printf("Copyright (C) 2013-2016  Tomas Racek (tom@krab1k.net)\n");
	printf("This program comes with ABSOLUTELY NO WARRANTY. This is free software, and you are welcome to redistribute it under certain conditions; see COPYING file for details.\n");
}

/* Prints help if -h/--help is issued */
static void print_help(void) {

	printf("Usage: neemp [OPTION] ...\n");
	printf("Perform EEM (Electronegativity Equalization Method) parameterization. Compute atomic charges using EEM.\n");

	printf("\nOptions:\n");
	printf("  -h, --help			 display this help and exit\n");
	printf("      --version			 display version information and exit\n");
	printf("      --max-threads N		 use up to N threads to solve EEM system in parallel\n");
	printf("  -m, --mode MODE		 set mode for the NEEMP. Valid choices are: info, params, charges, quality, cover (required)\n");
	printf("  -p, --params-method METHOD set optimization method used for calculation of parameters. Valid choices are: lr-full, lr-full-brent, de, gm (optional)\n");
	printf("      --sdf-file FILE		 SDF file (required)\n");
	printf("      --atom-types-by METHOD	 classify atoms according to the METHOD. Valid choices are: Element, ElemBond or User.\n");
	printf("      --list-omitted-molecules	 list names of molecules for which we don't have charges or parameters loaded (mode dependent).\n");
	printf("Options specific to mode: params using linear regression as calculation method\n");
	printf("      --chg-file FILE            FILE with ab-initio charges (required)\n");
	printf("      --chg-stats-out-file FILE  output charges statistics to the FILE\n");
	printf("      --random-seed VALUE        set random seed\n");
	printf("      --kappa-max MAX            set maximum value for kappa (required)\n");
	printf("      --kappa VALUE              use only one kappa VALUE for parameterization\n");
	printf("      --fs-precision VALUE       resolution for the full scan (required)\n");
	printf("      --kappa-preset PRESET      set kappa-max and fs-precision to safe values. Valid choices are: small, protein.\n");
	printf("Options specific to mode: params using optimization method (differential evolution, guided minimization)\n");
	printf("      --om-pop-size VALUE        set population size for optimization method (optional).\n");
	printf("      --om-iters COUNT  	     set the maximum number of iterations for optimization method (optional).\n");
	printf("      --om-threads      		 set number of threads for optimization method (optional).\n");
	printf("      --om-polish VALUE    		 apply local minimzation on parameters. Valid choices: 0 (off), 1 (result), 2 (during evolving), 3 (at the initial population). \n");
	printf("Options specific to mode: params using differential evolution\n");
	printf("      --de-f VALUE               set mutation constant for DE (optional).\n");
	printf("      --de-cr VALUE              set crossover recombination constant for DE (optional).\n");
	printf("      --de-dither                set the mutation constant to random value from [0.5;1] for ech iteration (optional).\n");
	printf("      --de-fix-kappa      		 set kappa to one fixed value (optional).\n");
	printf("Options specific to mode: params using guided minimization\n");
	printf("      --gm-iterations-beg  		 set number of minimization iterations for each reasonable vector of parameters (optional).\n");
	printf("      --gm-iterations-end  		 set number of minimization itertions for the best to polish the final result (optional).\n");
	printf("Other options:\n");
	printf("      --par-out-file FILE        output the parameters to the FILE\n");
	printf("  -d, --discard METHOD           perform discarding with METHOD. Valid choices are: iterative, simple and off. Default is off.\n");
	printf("  -s, --sort-by STAT             sort solutions by STAT. Valid choices are: R, R2, R_w, spearman, RMSD, RMSD_avg, D_max, D_avg. We strongly advise using R_w for method DE.\n");
	printf("      --limit-iters COUNT        set the maximum number of iterations for discarding.\n");
	printf("      --limit-time HH:MM:SS      set the maximum time for discarding in format hours:minutes:seconds.\n");
	printf("      --check-charges      	 warn about molecules with abnormal differences between QM and EEM charges.\n");
	printf("Options specific to mode: charges\n");
	printf("      --par-file FILE		 FILE with EEM parameters (required)\n");
	printf("      --chg-out-file FILE	 Output charges to the FILE (required)\n");

	printf("\nExamples:\n");
	printf("neemp -m info --sdf-file molecules.sdf --atom-types-by Element\n\
		Display information about the training set in the file molecules.sdf. Group atoms according to the elements only.\n");

	printf("neemp -m params --sdf-file molecules.sdf --chg-file charges.chg --kappa-max 1.0 --fs-precision 0.2 --sort-by RMSD.\n\
		Compute parameters for the given molecules in file molecules.sdf and ab-initio charges in charges.chg. Set maximum value for kappa to 1.0, step for the full scan to 0.2, no iterative refinement, sort results according to the relative mean square deviation.\n");
	printf("neemp -m params -p gm --sdf-file molecules.sdf --chg-file charges.chg --om-pop-size 50 -gm-iterations-beg 1000 -gm-iterations-end 500 --random-seed 1234 -vv.\n\
		Compute parameters for the given molecules in file molecules.sdf and ab-initio charges in charges.chg. The chosen optimization method: guided minimization will create 250 vectors (each vector consists of all parameters) and minimized reasonably good ones for 1000 iterations. The best of them will be minimized again, for 500 iterations.\n");

	printf("neemp -m charges --sdf-file molecules.sdf --par-file parameters --chg-out-file output.chg\n\
		Calculate and store EEM charges to the file output.chg\n");
}

/* Parse command line options */
void parse_options(int argc, char **argv) {

	int c;
	int option_idx;

	while((c = getopt_long(argc, argv, "p:vd:s:hm:", long_options, &option_idx)) != -1)
	{

		switch(c) {
			case 'h':
				print_help();
				exit(RETURN_OK);

			case 'm': /* mode */
				if(!strcmp(optarg, "info"))
					s.mode = MODE_INFO;
				else if (!strcmp(optarg, "charges"))
					s.mode = MODE_CHARGES;
				else if (!strcmp(optarg, "params"))
					s.mode = MODE_PARAMS;
				else if (!strcmp(optarg, "quality"))
					s.mode = MODE_QUALITY;
				else if (!strcmp(optarg, "cover"))
					s.mode = MODE_COVER;
				else
					EXIT_ERROR(ARG_ERROR, "Invalid mode: %s\n", optarg);
				break;

			case 'p': /* parameters' calculation optimization method */
				if (!strcmp(optarg, "lr-full"))
					s.params_method = PARAMS_LR_FULL;
				else if (!strcmp(optarg, "lr-full-brent"))
					s.params_method = PARAMS_LR_FULL_BRENT;
				else if (!strcmp(optarg, "de"))
					s.params_method = PARAMS_DE;
				else if (!strcmp(optarg, "gm"))
					s.params_method = PARAMS_GM;
				else 
					EXIT_ERROR(ARG_ERROR, "Invalid params-method: %s\n", optarg);
				break;

			case 'v':
				s.verbosity++;
				break;

			case 'd': /* discard */
				if(!strcmp(optarg, "off"))
					s.discard = DISCARD_OFF;
				else if(!strcmp(optarg, "iterative"))
					s.discard = DISCARD_ITER;
				else if(!strcmp(optarg, "simple"))
					s.discard = DISCARD_SIMPLE;
				else
					EXIT_ERROR(ARG_ERROR, "Invalid discarding mode: %s\n", optarg);
				break;

			case 's': /* sort-by */
				if(!strcmp(optarg, "R"))
					s.sort_by = SORT_R;
				else if(!strcmp(optarg, "R2"))
					s.sort_by = SORT_R2;
				else if(!strcmp(optarg, "R_w"))
					s.sort_by = SORT_RW;
				else if(!strcmp(optarg, "Spearman"))
					s.sort_by = SORT_SPEARMAN;
				else if(!strcmp(optarg, "RMSD"))
					s.sort_by = SORT_RMSD;
				else if(!strcmp(optarg, "RMSD_avg"))
					s.sort_by = SORT_RMSD_AVG;
				else if(!strcmp(optarg, "D_avg"))
					s.sort_by = SORT_D_AVG;
				else if(!strcmp(optarg, "D_max"))
					s.sort_by = SORT_D_MAX;
				else
					EXIT_ERROR(ARG_ERROR, "Invalid sort-by value: %s\n", optarg);
				break;

			case 129:
				print_version();
				exit(RETURN_OK);

			case 130:
				strncpy(s.sdf_file, optarg, MAX_PATH_LEN - 1);
				break;

			case 131:
				strncpy(s.par_file, optarg, MAX_PATH_LEN - 1);
				break;

			case 132:
				strncpy(s.chg_file, optarg, MAX_PATH_LEN - 1);
				break;

			case 133:
				strncpy(s.chg_out_file, optarg, MAX_PATH_LEN - 1);
				break;

			case 134:
				strncpy(s.chg_stats_out_file, optarg, MAX_PATH_LEN - 1);
				break;

			case 135:
				strncpy(s.par_out_file, optarg, MAX_PATH_LEN - 1);
				break;

			case 136:
				strncpy(s.atb_file, optarg, MAX_PATH_LEN - 1);
				break;

			case 137:
				s.random_seed = atoi(optarg);
				break;

			case 140:
				s.kappa_max = (float) atof(optarg);
				break;

			case 141:
				s.kappa_set = (float) atof(optarg);
				break;
			case 142:
				if(!strcmp(optarg, "small")) {
					s.kappa_max = 1.5f;
					s.full_scan_precision = 0.1f;
				}
				else if(!strcmp(optarg, "protein")) {
					s.kappa_max = 0.01f;
					s.full_scan_precision = 0.001f;
				}
				else
					EXIT_ERROR(ARG_ERROR, "Invalid kappa-preset value: %s\n", optarg);
				break;
			case 143:
				s.full_scan_precision = (float) atof(optarg);
				break;
			case 151: /* at-customization */
				if(!strcmp(optarg, atom_types_by_strings[AT_CUSTOM_ELEMENT]))
					s.at_customization = AT_CUSTOM_ELEMENT;
				else if(!strcmp(optarg, atom_types_by_strings[AT_CUSTOM_ELEMENT_BOND]))
					s.at_customization = AT_CUSTOM_ELEMENT_BOND;
				else if(!strcmp(optarg, atom_types_by_strings[AT_CUSTOM_USER]))
					s.at_customization = AT_CUSTOM_USER;
				else
					EXIT_ERROR(ARG_ERROR, "Invalid atom-type-by value: %s\n", optarg);
				break;

			case 152:
				s.tabu_size = (float) atof(optarg);
				break;
			case 160:
				s.limit_iters =  atoi(optarg);
				break;
			case 161: {
						 char *part;
						 int hours = 0;
						 int mins = 0;
						 int secs = 0;

						 part = strtok(optarg, ":");
						 if(part != NULL)
							 hours = atoi(part);
						 part = strtok(NULL, ":");
						 if(part != NULL)
							 mins = atoi(part);
						 part = strtok(NULL, ":");
						 if(part != NULL)
							 secs = atoi(part);

						 s.limit_time = 3600 * hours + 60 * mins + secs;
						 break;
					 }
			case 170:
					 s.check_charges = 1;
					 break;
			case 171:
					 s.max_threads =  atoi(optarg);
					 break;
			case 172:
					 s.list_omitted_molecules = 1;
					 break;
			case 173:
					 s.extra_precise = 1;
					 break;
			/* DE settings */
			case 180:
					 s.population_size = atoi(optarg);
					 break;
			case 181:
					 s.mutation_constant = (float) atof(optarg);
					 break;
			case 182:
					 s.recombination_constant = (float) atof(optarg);
					 break;
			case 183:
					 s.om_iters = atoi(optarg);
					 break;
			case 184: {
						 char *part;
						 int hours = 0;
						 int mins = 0;
						 int secs = 0;

						 part = strtok(optarg, ":");
						 if(part != NULL)
							 hours = atoi(part);
						 part = strtok(NULL, ":");
						 if(part != NULL)
							 mins = atoi(part);
						 part = strtok(NULL, ":");
						 if(part != NULL)
							 secs = atoi(part);

						 s.om_time = 3600 * hours + 60 * mins + secs;
						 break;
					 }
			case 186:
					 s.dither = 1;
					 break;
			case 188:
					 s.fixed_kappa = (float)atof(optarg);
					 break;
			case 189:
					 s.om_threads = atoi(optarg);
					 break;
			case 190:
					  s.polish = atoi(optarg);
					  break;
			/* GM settings */
			case 192:
					  s.gm_iterations_beg = atoi(optarg);
					  break;
			case 193:
					  s.gm_iterations_end = atoi(optarg);
					  break;
			case '?':
				EXIT_ERROR(ARG_ERROR, "%s", "Try -h/--help.\n");
			default:
				EXIT_ERROR(ARG_ERROR, "%s", "We should not be here!\n");
		}

	}
}

/* Check if options are set correctly */
void check_settings(void) {

	if(s.mode == MODE_NOT_SET)
		EXIT_ERROR(ARG_ERROR, "%s", "No mode set. Use '-m MODE', where MODE = info, params, quality, charges, cover.\n");

	if(s.sdf_file[0] == '\0')
		EXIT_ERROR(ARG_ERROR, "%s", "No .sdf file provided. Use '--sdf-file FILE'.\n");

	if(s.max_threads < 1)
		EXIT_ERROR(ARG_ERROR, "%s", "Maximum number of threads has to be at least 1 (default).\n");

	if(s.om_threads < 1)
		EXIT_ERROR(ARG_ERROR, "%s", "Maximum number of OM threads has to be at least 1 (default).\n");

	if(s.max_threads < s.om_threads)
		EXIT_ERROR(ARG_ERROR, "%s", "Maximum number of OM threads has to be smaller than maximum number of threads.\n");

	if(s.mode == MODE_PARAMS) {
		if(s.chg_file[0] == '\0')
			EXIT_ERROR(ARG_ERROR, "%s", "No .chg file provided. Use '--chg-file FILE'.\n");
		/* If user did not specify the optimization method for parameters calculation, set linear regression */
		if (s.params_method == PARAMS_NOT_SET)
			s.params_method = PARAMS_LR_FULL;

		if (s.params_method == PARAMS_LR_FULL || s.params_method == PARAMS_LR_FULL_BRENT) {
			if(s.full_scan_precision < 0)
				EXIT_ERROR(ARG_ERROR, "%s", "Full scan precision must greater than zero.\n");

			if(s.kappa_max < 0)
				EXIT_ERROR(ARG_ERROR, "%s", "Maximum for kappa must be greater than zero.\n");

			if(s.full_scan_precision > s.kappa_max)
				EXIT_ERROR(ARG_ERROR, "%s", "Full scan precision must be less than kappa max.\n");
		}
		
		if (s.params_method == PARAMS_DE) {
			/* All settings are optional, so check for mistakes and set defaults */
			if (s.population_size <1)
				s.population_size = 1000; /* 1.2 * (ts.atom_types_count * 2 + 1); */
			if (s.om_iters == NO_LIMIT_ITERS && s.om_time == NO_LIMIT_TIME)
				s.om_iters = 2000;
			if (s.mutation_constant < 0) /* If not set */
				s.mutation_constant = 0.75;
			if (s.recombination_constant < 0)
				s.recombination_constant = 0.7;
			if (s.polish == -1)
				s.polish = 3;
			if (s.sort_by == SORT_NOT_SET)
				s.sort_by = SORT_RMSD_AVG;
		}

		if (s.params_method == PARAMS_GM) {
			/* All settings are optional, so check for mistakes */
			if (s.population_size < 1)
				s.population_size = 100; /* 1.2 * (ts.atom_types_count * 2 + 1); */
			if (s.gm_iterations_beg < 1 || s.gm_iterations_end < 1)
				EXIT_ERROR(ARG_ERROR, "%s", "Number of minimization iterations for GM has to be positive.\n");
			if (s.sort_by == SORT_NOT_SET)
				s.sort_by = SORT_RMSD_AVG;
		}

		if (s.random_seed == -1)
			s.random_seed = 123;

		if(s.tabu_size < 0.0f || s.tabu_size > 1.0f)
			EXIT_ERROR(ARG_ERROR, "%s", "Tabu size has to be number in range [0.0; 1.0]\n");

		if(s.limit_iters != NO_LIMIT_ITERS && s.limit_iters > 100000)
			EXIT_ERROR(ARG_ERROR, "%s", "Number of iterations should be no higher than 1e6.\n");

		if(s.limit_time != NO_LIMIT_TIME && s.limit_time > 36000 * 1000)
			EXIT_ERROR(ARG_ERROR, "%s", "Maximum time should not be higher than 1000 hours.\n");

		if ((s.params_method == PARAMS_LR_FULL || s.params_method == PARAMS_LR_FULL_BRENT) && s.sort_by == SORT_NOT_SET)
			s.sort_by = SORT_R2;

		/* TODO verify with Tomas if this is the intended behavior */
		if(s.params_method == PARAMS_LR_FULL_BRENT /*!s.full_scan_only*/ && (s.sort_by != SORT_R && s.sort_by != SORT_R2 && s.sort_by != SORT_SPEARMAN && s.sort_by != SORT_RW))
			EXIT_ERROR(ARG_ERROR, "%s", "Full scan must be used for sort-by other than R, R2 or Spearman.\n");

	} else if(s.mode == MODE_CHARGES) {
		if(s.par_file[0] == '\0')
			EXIT_ERROR(ARG_ERROR, "%s", "No .par file provided. Use option '--par-file'.\n");

		if(s.chg_out_file[0] == '\0')
			EXIT_ERROR(ARG_ERROR, "%s", "No .chg output file provided. Use option '--chg-out-file'.\n");
	} else if(s.mode == MODE_QUALITY) {
		if(s.chg_file[0] == '\0')
			EXIT_ERROR(ARG_ERROR, "%s", "No .chg file provided. Use option '--chg-file'.\n");

		if(s.par_file[0] == '\0')
			EXIT_ERROR(ARG_ERROR, "%s", "No .par file provided. Use option '--par-file'.\n");
	} else if(s.mode == MODE_COVER) {
		if(s.par_file[0] == '\0')
			EXIT_ERROR(ARG_ERROR, "%s", "No .par file provided. Use option '--par-file'.\n");
	}

	if(s.at_customization == AT_CUSTOM_USER && s.atb_file[0] == '\0')
		EXIT_ERROR(ARG_ERROR, "%s", "File with user defined types (option '--atb-file') must be provided when runned with '--atom-types-by User'\n");
}

void print_settings(void) {

	printf("\nSettings:\n\n");

	printf("Mode: ");
	switch(s.mode) {
		case MODE_INFO:
			printf("info (print info about the training set)\n");
			break;
		case MODE_PARAMS:
			printf("params (calculate EEM parameters)");
			if (s.params_method == PARAMS_LR_FULL)
				printf(" with full scan of kappa range with linear regression\n");
			if (s.params_method == PARAMS_LR_FULL_BRENT)
				printf(" with full scan of kappa with Brent's method\n");
			if (s.params_method == PARAMS_DE)
				printf(" with differential evolution method\n");
			break;
		case MODE_CHARGES:
			printf("charges (calculate EEM charges)\n");
			break;
		case MODE_QUALITY:
			printf("quality (perform quality validation of EEM parameters)\n");
			break;
		case MODE_COVER:
			printf("cover (perform coverage validation of EEM parameters)\n");
			break;
		case MODE_NOT_SET:
			assert(0);
	}

	printf("\nFiles:\n");
	printf(" Structural (.sdf) file: %s\n", s.sdf_file);
	if(s.par_file[0] != '\0')
		printf(" Parameters (.par) file: %s\n", s.par_file);

	if(s.chg_file[0] != '\0')
		printf(" Charges (.chg) file: %s\n", s.chg_file);

	if(s.atb_file[0] != '\0')
		printf(" Atom types (.atb) file: %s\n", s.atb_file);

	if(s.par_out_file[0] != '\0')
		printf(" Parameters (.par) output file: %s\n", s.par_out_file);

	if(s.chg_out_file[0] != '\0')
		printf(" Charges output (.chg) file: %s\n", s.chg_out_file);

	if(s.chg_stats_out_file[0] != '\0')
		printf(" Charges stats output (.chgs) file: %s\n", s.chg_stats_out_file);

	printf("\nAtom types grouped by: ");
	switch(s.at_customization) {
		case AT_CUSTOM_ELEMENT:
			printf("%s (element)\n", atom_types_by_strings[AT_CUSTOM_ELEMENT]);
			break;
		case AT_CUSTOM_ELEMENT_BOND:
			printf("%s (element + bond order)\n", atom_types_by_strings[AT_CUSTOM_ELEMENT_BOND]);
			break;
		case AT_CUSTOM_USER:
			printf("%s (from external file)\n", atom_types_by_strings[AT_CUSTOM_USER]);
			break;
	}
    printf("\nMaximum number of threads: %d\n", s.max_threads);
	if (s.mode == MODE_PARAMS && s.params_method == PARAMS_DE)
		printf("Maximum number of threads used for DE: %d\n", s.om_threads);
	printf("\nVerbosity level: ");
	switch(s.verbosity) {
		case VERBOSE_MINIMAL:
			printf("0 (only minimal info)\n");
			break;
		case VERBOSE_DISCARD:
			printf("1 (include discarding info)\n");
			break;
		case VERBOSE_KAPPA:
			printf("2 (include discarding + kappa search/de info)\n");
			break;
		default:
			printf("%u (I won't tell you more than on level 2. Sorry about that.)\n", s.verbosity);
			break;
	}

	if(s.mode == MODE_PARAMS) {
		printf("\nSort by: ");
		switch(s.sort_by) {
			case SORT_R:
				printf("R");
				break;
			case SORT_R2:
				printf("R2");
				break;
			case SORT_RW:
				printf("RW");
				break;
			case SORT_SPEARMAN:
				printf("Sp");
				break;
			case SORT_RMSD:
				printf("RMSD");
				break;
			case SORT_RMSD_AVG:
				printf("average RMSD");
				break;
			case SORT_D_AVG:
				printf("Avg D");
				break;
			case SORT_D_MAX:
				printf("Max D");
				break;
			case SORT_NOT_SET:
				printf("not set");
				break;
		}

		if(s.sort_by == SORT_R || s.sort_by == SORT_R2 || s.sort_by == SORT_RW || s.sort_by == SORT_SPEARMAN)
			printf(" (higher is better)\n");
		else
			printf(" (lower is better)\n");

		printf("\nDiscarding:\n");
		printf(" Mode: ");
		switch(s.discard) {
			case DISCARD_OFF:
				printf("off\n");
				break;
			case DISCARD_SIMPLE:
				printf("simple (try each molecule once)\n");
				break;
			case DISCARD_ITER:
				printf("iterative (try molecules at random until limit is reached)\n");
				break;
		}
		printf(" Time limit: ");
		if(s.limit_time == NO_LIMIT_TIME)
			printf("none set\n");
		else {

			unsigned tmp = (unsigned) s.limit_time;
			unsigned secs = tmp % 60;
			tmp /= 60;
			unsigned mins = tmp % 60;
			tmp /= 60;
			unsigned hours = tmp;
			printf("%u hour(s) %u min(s) %u sec(s)\n", hours, mins, secs);


		}
		printf(" Iterations limit: ");
		if(s.limit_iters == NO_LIMIT_ITERS)
			printf("none set\n");
		else
			printf("%d\n", s.limit_iters);

		if (s.params_method == PARAMS_LR_FULL || s.params_method == PARAMS_LR_FULL_BRENT) {
			printf("\nKappa search:\n");
			printf(" Mode: ");
			if(s.kappa_set > 0.0f)
				printf("single kappa value = %5.3f\n", s.kappa_set);
			else {
				printf("full scan from 0.0 to %5.3f with step %5.3f", s.kappa_max, s.full_scan_precision);

				if(s.params_method == PARAMS_LR_FULL_BRENT)
					printf(" + iterative refinement\n");
				else
					printf("\n");
			}
		}
		if (s.params_method == PARAMS_DE) {
			printf("\n Differential evolution settings:\n");
			printf("\t - population size %d\n", s.population_size);
			printf("\t - max iterations  %d\n", s.om_iters);
			if (s.polish != 0) {
				printf("\t - polishing ");
				if (s.polish >= 1)
					printf("the result");
				if (s.polish >= 2)
					printf(" and during evolving");
				if (s.polish >= 3)
					printf(" and some of the initial population");
				printf("\n");
			}
			printf("\t - mutation constant %5.3lf\n", s.mutation_constant);
			printf("\t - recombination constant %5.3lf\n", s.recombination_constant);
			if (s.dither != 0)
				printf("\t - dither on\n");
			if (s.fixed_kappa > 0)
				printf("\t - kappa fixed on value %5.3lf\n", s.fixed_kappa);


		}

		if (s.params_method == PARAMS_GM) {
			printf("\nGuided minimization settings:\n");
			printf("\t - set size %d\n", s.population_size);
			printf("\t - iterations for set at the beginning %d\n", s.gm_iterations_beg);
			printf("\t - iterations for the result at the end %d\n", s.gm_iterations_end);
			printf("\t - threads used for minimization %d\n", s.om_threads);
		}
	}

	printf("\n");
}
