/*
 * NEEMP - tabu.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013, 2014
 *
 * */

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
