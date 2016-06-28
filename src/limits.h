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
