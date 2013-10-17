/*
 * NEEMP - settings.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#include <getopt.h>
#include <string.h>

#include "neemp.h"
#include "settings.h"

static void print_help(void);

extern struct settings s;

static struct option long_options[] = {
	{"help", no_argument, 0, 'h'},
	{"mode", required_argument, 0, 'm'},
	{"fs-only", no_argument, 0, 'f'},
	{"sdf-file", required_argument, 0, 10},
	{"par-file", required_argument, 0, 11},
	{"chg-file", required_argument, 0, 12},
	{"chg-out-file", required_argument, 0, 13},
	{"chg-stats-out-file", required_argument, 0, 14},
	{"par-out-file", required_argument, 0, 15},
	{"kappa-max", required_argument, 0, 20},
	{"kappa", required_argument, 0, 21},
	{"fs-precision", required_argument, 0, 22},
	{"sort-by", required_argument, 0, 30},
	{"at-customization", required_argument, 0, 31},
	{NULL, 0, 0, 0}
};

/* Initialize default settings */
void s_init(void) {

	memset(s.sdf_file, 0x0, MAX_PATH_LEN * sizeof(char));
	memset(s.chg_file, 0x0, MAX_PATH_LEN * sizeof(char));
	memset(s.par_file, 0x0, MAX_PATH_LEN * sizeof(char));
	memset(s.par_out_file, 0x0, MAX_PATH_LEN * sizeof(char));
	memset(s.chg_out_file, 0x0, MAX_PATH_LEN * sizeof(char));
	memset(s.chg_stats_out_file, 0x0, MAX_PATH_LEN * sizeof(char));

	s.mode = MODE_NOT_SET;
	s.full_scan_only = 0;
	s.full_scan_precision = 0.0f;
	s.kappa_max = 0.0f;
	s.full_scan_precision = 0.0f;
	s.sort_by = SORT_R;
	s.at_customization = AT_CUSTOM_ELEMENT_BOND;
}

/* Prints help if -h/--help is issued */
static void print_help(void) {

	printf("%s %s\n\n", APP_NAME, APP_VERSION);
}

/* Parse command line options */
void parse_options(int argc, char **argv) {

	int c;
	int option_idx;

	while(1) {
		c = getopt_long(argc, argv, "fhm:", long_options, &option_idx);
		if(c == -1)
			break;

		switch(c) {
			case 'f':
				s.full_scan_only = 1;
				break;
			case 'h':
				print_help();
				exit(RETURN_OK);
			case 'm':
				if(!strcmp(optarg, "info"))
					s.mode = MODE_INFO;
				else if (!strcmp(optarg, "charges"))
					s.mode = MODE_CHARGES;
				else if (!strcmp(optarg, "params"))
					s.mode = MODE_PARAMS;
				else
					EXIT_ERROR(ARG_ERROR, "Invalid mode: %s\n", optarg);
				break;
			case 10:
				strncpy(s.sdf_file, optarg, MAX_PATH_LEN - 1);
				break;
			case 11:
				strncpy(s.par_file, optarg, MAX_PATH_LEN - 1);
				break;
			case 12:
				strncpy(s.chg_file, optarg, MAX_PATH_LEN - 1);
				break;
			case 13:
				strncpy(s.chg_out_file, optarg, MAX_PATH_LEN - 1);
				break;
			case 14:
				strncpy(s.chg_stats_out_file, optarg, MAX_PATH_LEN - 1);
				break;
			case 15:
				strncpy(s.par_out_file, optarg, MAX_PATH_LEN - 1);
				break;
			case 20:
				s.kappa_max = (float) atof(optarg);
				break;
			case 21:
				s.kappa_set = (float) atof(optarg);
				break;
			case 22:
				s.full_scan_precision = (float) atof(optarg);
				break;
			case 30:
				if(!strcmp(optarg, "R"))
					s.sort_by = SORT_R;
				else if(!strcmp(optarg, "RMSD"))
					s.sort_by = SORT_RMSD;
				else if(!strcmp(optarg, "MSE"))
					s.sort_by = SORT_MSE;
				else if(!strcmp(optarg, "D_avg"))
					s.sort_by = SORT_D_AVG;
				else if(!strcmp(optarg, "D_max"))
					s.sort_by = SORT_D_MAX;
				else
					EXIT_ERROR(ARG_ERROR, "Invalid sort-by value: %s\n", optarg);
				break;
			case 31:
				if(!strcmp(optarg, "element"))
					s.at_customization = AT_CUSTOM_ELEMENT;
				else if(!strcmp(optarg, "element_bond"))
					s.at_customization = AT_CUSTOM_ELEMENT_BOND;
				else if(!strcmp(optarg, "partner"))
					s.at_customization = AT_CUSTOM_PARTNER;
				else if(!strcmp(optarg, "valence"))
					s.at_customization = AT_CUSTOM_VALENCE;
				else
					EXIT_ERROR(ARG_ERROR, "Invalid at-customization value: %s\n", optarg);
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

	if(s.at_customization == AT_CUSTOM_PARTNER || s.at_customization == AT_CUSTOM_VALENCE)
		EXIT_ERROR(ARG_ERROR, "%s", "These atom type customizations are not implemented right now. Stay tuned.\n");

	if(s.mode == MODE_NOT_SET)
		EXIT_ERROR(ARG_ERROR, "%s", "No mode set.\n");

	if(s.mode == MODE_PARAMS) {
		if(s.sdf_file[0] == '\0')
			EXIT_ERROR(ARG_ERROR, "%s", "No .sdf file provided.\n");
		if(s.chg_file[0] == '\0')
			EXIT_ERROR(ARG_ERROR, "%s", "No .chg file provided.\n");
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
			if(s.full_scan_only)
				EXIT_ERROR(ARG_ERROR, "%s", "Cannot use full scan if single kappa is selected.\n");
		}
	} else if(s.mode == MODE_CHARGES) {
		if(s.sdf_file[0] == '\0')
			EXIT_ERROR(ARG_ERROR, "%s", "No .sdf file provided.\n");
		if(s.par_file[0] == '\0')
			EXIT_ERROR(ARG_ERROR, "%s", "No .par file provided.\n");
		if(s.chg_out_file[0] == '\0')
			EXIT_ERROR(ARG_ERROR, "%s", "No .chg output file provided.\n");
	}
}
