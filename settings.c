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
	{"sdf-file", required_argument, 0, 10},
	{"par-file", required_argument, 0, 11},
	{"chg-file", required_argument, 0, 12}
};


/* Prints help if -h/--help is issued */
static void print_help(void) {

	printf("%s %s\n\n", APP_NAME, APP_VERSION);
}

/* Parse command line options */
void parse_options(int argc, char **argv) {

	int c;
	int option_idx;

	while(1) {
		c = getopt_long(argc, argv, "h", long_options, &option_idx);
		if(c == -1)
			break;

		switch(c) {
			case 'h':
				print_help();
				exit(RETURN_OK);
			case 10:
				strncpy(s.sdf_filename, optarg, MAX_PATH_LEN);
				s.sdf_filename[MAX_PATH_LEN - 1] = '\0';
				break;
			case 11:
				strncpy(s.par_filename, optarg, MAX_PATH_LEN);
				s.par_filename[MAX_PATH_LEN - 1] = '\0';
				break;
			case 12:
				strncpy(s.chg_filename, optarg, MAX_PATH_LEN);
				s.chg_filename[MAX_PATH_LEN - 1] = '\0';
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

}
