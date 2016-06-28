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
