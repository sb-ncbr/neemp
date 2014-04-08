/*
 * NEEMP - config.h
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013, 2014
 *
 * */

#ifndef __CONFIG_H__
#define __CONFIG_H__

/* Constraints */
#define MAX_PATH_LEN 100
#define MAX_LINE_LEN 80

#define MAX_MOLECULES 1E5
#define MAX_ATOM_TYPES 100

#define MIN_ATOMS_PER_MOLECULE 2
#define MAX_ATOMS_PER_MOLECULE 1E5
#define MIN_BONDS_PER_MOLECULE 1
#define MAX_BONDS_PER_MOLECULE 1E5

/* Default per molecule warning constants for --check-charges */
#define WARN_MIN_R 0.2F
#define WARN_MAX_RMSD 0.5F
#define WARN_MAX_MSE 1.0F
#define WARN_MAX_D_AVG 0.2F
#define WARN_MAX_D_MAX 1.0F

#define WARN_MIN_RCOND 0.000005F

#endif /* __CONFIG_H__ */
