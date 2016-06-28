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
#include <stdlib.h>

#include "neemp.h"
#include "tabu.h"

void t_init(struct tabu * const t, int size) {

	assert(t != NULL);

	t->size = size;

	t->list = (int *) malloc(sizeof(int) * size);
	if(!t->list)
		EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for tabu list.\n");

	for(int i = 0; i < t->size; i++)
		t->list[i] = -1;

	t->current_idx = 0;
}

int t_is_banned(const struct tabu * const t, int idx) {

	assert(t != NULL);

	for(int i = 0; i < t->size; i++)
		if(t->list[i] == idx)
			return 1;

	return 0;
}

void t_update(struct tabu * const t, int idx) {

	assert(t != NULL);

	if(t->size != 0) {
		t->list[t->current_idx] = idx;
		t->current_idx = (t->current_idx + 1) % t->size;
	}
}

void t_destroy(struct tabu * const t) {

	assert(t != NULL);

	free(t->list);
}
