/*
 * NEEMP - limits.h
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#ifndef __LIMITS_H__
#define __LIMITS_H__

#include <time.h>

#define NO_LIMIT_ITERS -1
#define NO_LIMIT_TIME -1

struct limit {

	int iters_max;
	int iters_current;

	time_t time_max;
};

void l_init(struct limit *l, int limit_iters, time_t limit_time);
int l_check(struct limit *l);

#endif /* __LIMITS_H__ */
