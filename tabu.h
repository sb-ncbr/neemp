/*
 * NEEMP - tabu.h
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

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
