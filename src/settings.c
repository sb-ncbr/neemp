/*
 * NEEMP - settings.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013, 2014
 *
 * */

#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <string.h>

#include "limits.h"
#include "neemp.h"
#include "settings.h"
#include "latin_random.h"

extern struct settings s;

static void print_help(void);
static void print_version(void);

static char *atom_types_by_strings[] = {"Element", "ElemBond", "Partner", "Valence"};

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
	{"wgh-file", required_argument, 0, 136},
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
	{"rw", required_argument, 0, 173},
	{"de-pop-size", required_argument, 0, 180},
	{"de-f", required_argument, 0, 181},
	{"de-cr", required_argument, 0, 182},
	{"de-iters-max", required_argument, 0, 183},
	{"de-time-max", required_argument, 0, 184},
	{"de-take-only-best", no_argument, 0, 185},
	{"de-dither", no_argument, 0, 186},
	{"de-evolve-partially", no_argument, 0, 187},
	{"de-fix-kappa",required_argument, 0, 188},
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
	memset(s.par_out_file, 0x0, MAX_PATH_LEN * sizeof(char));
	memset(s.chg_out_file, 0x0, MAX_PATH_LEN * sizeof(char));
	memset(s.chg_stats_out_file, 0x0, MAX_PATH_LEN * sizeof(char));

	s.random_seed = -1;
	s.mode = MODE_NOT_SET;
	s.params_method = PARAMS_NOT_SET;
	s.full_scan_precision = 0.0f;
	s.kappa_max = 0.0f;
	s.full_scan_precision = 0.0f;
	s.population_size = 0;
	s.recombination_constant = 0.0f;
	s.mutation_constant = 0.0f;
	s.dither = 0;
	s.evolve_by_element = 0;
	s.fixed_kappa = 0;
	s.take_only_best = 0;
	s.limit_de_iters = NO_LIMIT_ITERS;
	s.limit_de_time = NO_LIMIT_TIME;
	s.sort_by = SORT_R2;
	s.at_customization = AT_CUSTOM_ELEMENT_BOND;
	s.discard = DISCARD_OFF;
	s.tabu_size = 0.0f;
	s.limit_iters = NO_LIMIT_ITERS;
	s.limit_time = NO_LIMIT_TIME;
	s.check_charges = 0;
	s.max_threads = 1;
	s.list_omitted_molecules = 0;
	s.rw = 0.5f;
}

/* Prints help if --version is issued */
static void print_version(void) {

	printf("%s %s\n", APP_NAME, APP_VERSION);
	printf("Copyright (C) 2013, 2014 Tomas Racek (tom@krab1k.net)\n");
	printf("Licence MIT: http://opensource.org/licenses/MIT\n");
	printf("This is free software: you are free to change and redistribute it.\n");
	printf("There is NO WARRANTY, to the extent permitted by law.\n");
}

/* Prints help if -h/--help is issued */
static void print_help(void) {

	printf("Usage: neemp [OPTION] ...\n");
	printf("Perform EEM (Electronegativity Equalization Method) parameterization. Compute atomic charges using EEM.\n");

	printf("\nOptions:\n");
	printf("  -h, --help			 display this help and exit\n");
	printf("      --version			 display version information and exit\n");
	printf("      --max-threads N		 use up to N threads to solve EEM system in parallel\n");
	printf("  -m, --mode MODE		 set mode for the NEEMP. Valid choices are: info, params, charges, cross, cover (required)\n");
	printf("  -p, --params-method METHOD set optimization method used for calculation of parameters. Valid choices are: lr-full, lr-full-brent, de (optional)\n");
	printf("      --sdf-file FILE		 SDF file (required)\n");
	printf("      --atom-types-by METHOD	 classify atoms according to the METHOD. Valid choices are: Element, ElemBond.\n");
	printf("      --list-omitted-molecules	 list names of molecules for which we don't have charges or parameters loaded (mode dependent).\n");
	printf("Options specific to mode: params using linear regression as calculation method\n");
	printf("      --chg-file FILE            FILE with ab-initio charges (required)\n");
	printf("      --chg-stats-out-file FILE  output charges statistics to the FILE\n");
	printf("      --random-seed VALUE        set random seed\n");
	printf("      --kappa-max MAX            set maximum value for kappa (required)\n");
	printf("      --kappa VALUE              use only one kappa VALUE for parameterization\n");
	printf("      --fs-precision VALUE       resolution for the full scan (required)\n");
	printf("      --kappa-preset PRESET      set kappa-max and fs-precision to safe values. Valid choices are: small, protein.\n");
	printf("Options specific to mode: params using differential evolution as calculation method\n");
	printf("      --de-pop-size VALUE        set population size for DE (optional).\n");
	printf("      --de-f VALUE               set mutation constant for DE (optional).\n");
	printf("      --de-cr VALUE              set crossover recombination constant for DE (optional).\n");
	printf("      --de-iters-max COUNT       set the maximum number of iterations for DE (optional).\n");
	printf("      --de-time-max HH:MM:SS     set the maximum time for DE in format hours:minutes:seconds (optional).\n");
	printf("      --de-take-only-best        turn on using only the better half of population for evolution.\n");
	printf("      --de-dither                set the mutation constant to random value from [0.5;1] for ech iteration, can improve convergence.\n");
	printf("      --de-evolve-partially      turn on evolution driven by sort per atom type.\n");
	printf("      --de-fix-kappa      		 set kappa to one fixed value.\n");
	printf("      --par-out-file FILE        output the parameters to the FILE\n");
	printf("  -d, --discard METHOD           perform discarding with METHOD. Valid choices are: iterative, simple and off. Default is off.\n");
	printf("  -s, --sort-by STAT             sort solutions by STAT. Valid choices are: R, R2, spearman, RMSD, D_max, D_avg.\n");
	printf("      --limit-iters COUNT        set the maximum number of iterations for discarding.\n");
	printf("      --limit-time HH:MM:SS      set the maximum time for discarding in format hours:minutes:seconds.\n");
	printf("      --check-charges      	 warn about molecules with abnormal differences between QM and EEM charges.\n");
	printf("Options specific to mode: charges\n");
	printf("      --par-file FILE		 FILE with EEM parameters (required)\n");
	printf("      --chg-out-file FILE	 Output charges to the FILE (required)\n");

	printf("\nExamples:\n");
	printf("neemp -m info --sdf-file molecules.sdf --atom-types-by Element\n\
		Display information about the training set in the file molecules.sdf. Group atoms according to the elements only.\n");

	printf("neemp -m params --sdf-file molecules.sdf --chg-file charges.chg --kappa-max 1.0 --fs-precision 0.2 --sort-by RMSD --fs-only.\n\
		Compute parameters for the given molecules in file molecules.sdf and ab-initio charges in charges.chg. Set maximum value for kappa to 1.0, step for the full scan to 0.2, no iterative refinement, sort results according to the relative mean square deviation.\n");
	printf("neemp -m params -p de --sdf-file molecules.sdf --chg-file charges.chg --sort-by R --de-pop-size 20 --de-iters-max 500 --de-evolve-partially -vv.\n\
		Compute parameters for the given molecules in file molecules.sdf and ab-initio charges in charges.chg. The chosen optimization method: differential evolution will create population of 20 sets of parameters and evolve these in maximum of 500 iterations. The fitness function evaluating the set of parameters is Pearson coefficient. Partial great improvements in evolution are permitted at the cost of slight decrease in total R.\n");

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
				else if (!strcmp(optarg, "cross"))
					s.mode = MODE_CROSS;
				else if (!strcmp(optarg, "cover"))
					s.mode = MODE_COVER;
				else
					EXIT_ERROR(ARG_ERROR, "Invalid mode: %s\n", optarg);
				break;

			case 'p': /*parameters' calculation optimization method */
				if (!strcmp(optarg, "lr-full"))
					s.params_method = PARAMS_LR_FULL;
				else if (!strcmp(optarg, "lr-full-brent"))
					s.params_method = PARAMS_LR_FULL_BRENT;
				else if (!strcmp(optarg, "de"))
					s.params_method = PARAMS_DE;
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
				strncpy(s.wgh_file, optarg, MAX_PATH_LEN - 1);
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
					 s.rw = (float) atof(optarg);
					 break;
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
					 s.limit_de_iters = atoi(optarg);
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

						 s.limit_de_time = 3600 * hours + 60 * mins + secs;
						 break;
					 }
			case 185:
					 s.take_only_best = 1;
					 break;
			case 186:
					 s.dither = 1;
					 break;
			case 187:
					 s.evolve_by_element = 1;
					 break;  
			case 188:
					 s.fixed_kappa = atof(optarg);
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
		EXIT_ERROR(ARG_ERROR, "%s", "No mode set.\n");

	if(s.sdf_file[0] == '\0')
		EXIT_ERROR(ARG_ERROR, "%s", "No .sdf file provided.\n");

	if(s.max_threads < 1)
		EXIT_ERROR(ARG_ERROR, "%s", "Maximum number of threads has to be at least 1 (default).\n");

	if(s.mode == MODE_PARAMS) {
		if(s.chg_file[0] == '\0')
			EXIT_ERROR(ARG_ERROR, "%s", "No .chg file provided.\n");
		//if user did not specify the optimization method for parameters calculation, set linear regression 
		if (s.params_method == PARAMS_NOT_SET)
			s.params_method = PARAMS_LR_FULL;

		if (s.params_method == PARAMS_LR_FULL || s.params_method == PARAMS_LR_FULL_BRENT) {
			if(s.kappa_set < 1e-10) {
				if(s.full_scan_precision < 1e-10)
					EXIT_ERROR(ARG_ERROR, "%s", "Full scan precision must be set correctly in params mode.\n");

				if(s.kappa_max < 1e-10)
					EXIT_ERROR(ARG_ERROR, "%s", "Maximum for kappa must be set correctly in params mode.\n");

				if(s.full_scan_precision > s.kappa_max)
					EXIT_ERROR(ARG_ERROR, "%s", "Full scan precision must be less than kappa max.\n");
			}
			else {
				if(s.full_scan_precision > 1e-10 || s.kappa_max > 1e-10)
					EXIT_ERROR(ARG_ERROR, "%s", "Cannot set full scan precision and/or kappa max if --kappa is used.\n");

				if(s.params_method == PARAMS_LR_FULL)
					EXIT_ERROR(ARG_ERROR, "%s", "Cannot use full scan if single kappa is selected.\n");
			}
		}
		
		if (s.params_method == PARAMS_DE) {
			if (s.mutation_constant == 0)
				s.mutation_constant = 0.75;
			if (s.recombination_constant == 0)
				s.recombination_constant = 0.7;
			if (s.population_size == 0)
				s.population_size = 20; //1.2*(ts.atom_types_count*2+1);
			if (s.limit_de_iters == NO_LIMIT_ITERS && s.limit_de_time == NO_LIMIT_TIME)
				s.limit_de_iters = 250;

		}

		if (s.random_seed == -1)
			s.random_seed = get_seed();

		if(s.tabu_size < 0.0f || s.tabu_size > 1.0f)
			EXIT_ERROR(ARG_ERROR, "%s", "Tabu size has to be number in range [0.0; 1.0]\n");

		if(s.limit_iters != NO_LIMIT_ITERS && s.limit_iters > 100000)
			EXIT_ERROR(ARG_ERROR, "%s", "Number of iterations should be no higher than 1e6.\n");

		if(s.limit_time != NO_LIMIT_TIME && s.limit_time > 36000 * 1000)
			EXIT_ERROR(ARG_ERROR, "%s", "Maximum time should not be higher than 1000 hours.\n");
		//TODO verify with Tomas if this is the intended behavior
		if(s.params_method == PARAMS_LR_FULL_BRENT /*!s.full_scan_only*/ && (s.sort_by != SORT_R && s.sort_by != SORT_R2 && s.sort_by != SORT_SPEARMAN && s.sort_by != SORT_RW))
			EXIT_ERROR(ARG_ERROR, "%s", "Full scan must be used for sort-by other than R, R2 or Spearman.\n");

	} else if(s.mode == MODE_CHARGES) {
		if(s.par_file[0] == '\0')
			EXIT_ERROR(ARG_ERROR, "%s", "No .par file provided.\n");

		if(s.chg_out_file[0] == '\0')
			EXIT_ERROR(ARG_ERROR, "%s", "No .chg output file provided.\n");
	} else if(s.mode == MODE_CROSS) {
		if(s.chg_file[0] == '\0')
			EXIT_ERROR(ARG_ERROR, "%s", "No .chg file provided.\n");

		if(s.par_file[0] == '\0')
			EXIT_ERROR(ARG_ERROR, "%s", "No .par file provided.\n");
	} else if(s.mode == MODE_COVER) {
		if(s.par_file[0] == '\0')
			EXIT_ERROR(ARG_ERROR, "%s", "No .par file provided.\n");
	}

	if(s.rw < 1.0f / 2.718281828f || s.rw >= 1.0f)
		EXIT_ERROR(ARG_ERROR, "%s", "--rw argument has to be in range [1/e; 1)\n");
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
		case MODE_CROSS:
			printf("cross (perform cross-validation of the EEM parameters)\n");
			break;
		case MODE_COVER:
			printf("cover (calculate parameters coverage)\n");
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

	if(s.wgh_file[0] != '\0')
		printf(" Weights (.wgh) file: %s\n", s.wgh_file);

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
	}

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
			case SORT_D_AVG:
				printf("Avg D");
				break;
			case SORT_D_MAX:
				printf("Max D");
				break;
		}

		if(s.sort_by == SORT_R || s.sort_by == SORT_R2 || s.sort_by == SORT_SPEARMAN)
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
	}

	printf("\n");
}
