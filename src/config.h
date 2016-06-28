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

#ifndef __CONFIG_H__
#define __CONFIG_H__

/* Constraints */
#define MAX_PATH_LEN 512
#define MAX_LINE_LEN 80

#define MAX_MOLECULES 50000000
#define MAX_ATOM_TYPES 300

#define MAX_ATOMS_PER_MOLECULE 1000000
#define MAX_BONDS_PER_MOLECULE 1000000

/* Default per molecule warning constants for --check-charges */
#define WARN_MIN_R 0.2F
#define WARN_MAX_RMSD 0.5F
#define WARN_MAX_MSE 1.0F
#define WARN_MAX_D_AVG 0.2F
#define WARN_MAX_D_MAX 1.0F
#define WARN_MAX_COND 50000

#define WARN_MIN_RCOND 0.000005F

#endif /* __CONFIG_H__ */
