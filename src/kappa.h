/*
 * NEEMP - kappa.h
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013, 2014
 *
 * */

#ifndef __KAPPA_H__
#define __KAPPA_H__

#include "subset.h"

void find_the_best_parameters_for_subset(struct subset * const ss);
void set_the_best(struct subset * const ss);

#endif /* __KAPPA_H__ */
