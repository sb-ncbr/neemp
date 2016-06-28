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

#ifndef __NEEMP_H__
#define __NEEMP_H__

#include <stdio.h>
#include <stdlib.h>

#define APP_NAME "NEEMP"
#define APP_VERSION "2.0"

#define RETURN_OK 0
#define ARG_ERROR 1
#define MEM_ERROR 2
#define IO_ERROR  3
#define RUN_ERROR 4

#define NOT_FOUND -1

#define EPS 10E-5

#define EXIT_ERROR(ERROR_CODE, fmt, ...) do { \
fprintf(stderr, "Error: " fmt, __VA_ARGS__);  \
if(ERROR_CODE == ARG_ERROR)\
fprintf(stderr, "Run 'neemp --help' to get a list of options.\n");\
exit(ERROR_CODE);                             \
} while(0)

#endif /* __NEEMP_H__ */
