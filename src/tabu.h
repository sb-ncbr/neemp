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

#ifndef __TABU_H__
#define __TABU_H__

struct tabu {

	int *list;
	int size;
	int current_idx;
};

void t_init(struct tabu * const t, int size);
void t_update(struct tabu * const t, int idx);
void t_destroy(struct tabu * const t);
int t_is_banned(const struct tabu * const t, int idx);

#endif /* __TABU_H__ */
