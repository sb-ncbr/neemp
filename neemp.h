/*
 * NEEMP - neemp.h
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#ifndef __NEEMP_H__
#define __NEEMP_H__

#include <stdio.h>
#include <stdlib.h>

#define APP_NAME "NEEMP"
#define APP_VERSION "1.1-git"

#define RETURN_OK 0
#define MEM_ERROR 1
#define IO_ERROR  2
#define RUN_ERROR 3

#define NOT_FOUND -1

#define EXIT_ERROR(ERROR_CODE, fmt, ...) do { \
fprintf(stderr, "Error: " fmt, __VA_ARGS__);  \
exit(ERROR_CODE);                             \
} while(0)

#endif /* __NEEMP_H__ */
