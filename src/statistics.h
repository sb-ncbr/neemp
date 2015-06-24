/*
 * NEEMP - statistics.h
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013, 2014
 *
 * */

#ifndef __STATISTICS_H__
#define __STATISTICS_H__

#include "subset.h"

void adjust_weights(struct subset * const ss);
void calculate_statistics(struct subset * const ss, struct kappa_data * const kd);
void check_charges(const struct kappa_data * const kd);

#endif /* __STATISTICS_H__ */
