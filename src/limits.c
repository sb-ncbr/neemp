/*
 * NEEMP - limits.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013, 2014
 *
 * */

#include <assert.h>
#include <time.h>

#include "limits.h"

/* Initialize the limit structure */
void l_init(struct limit *l, int limit_iters, time_t limit_time) {

	assert(l != NULL);

	l->iters_max = limit_iters;
	l->iters_current = 0;

	l->time_max = limit_time == NO_LIMIT_TIME ? NO_LIMIT_TIME : time(NULL) + limit_time;
}

/* Check if limits were reached */
int l_check(struct limit *l) {

	assert(l != NULL);

	if(l->iters_max != NO_LIMIT_ITERS && l->iters_current >= l->iters_max)
		return 1;

	if(l->time_max != NO_LIMIT_TIME && time(NULL) >= l->time_max)
		return 1;

	return 0;
}
