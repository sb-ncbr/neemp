/*
 * NEEMP - guidedmin.h
 *
 * by Jana Pazurikova (pazurikova@ics.muni.cz)
 * 2016
 *
 * */

#ifndef __GUIDEDMIN_H__
#define __GUIDEDMIN_H__

#include "subset.h"
#include "diffevolution.h"

void run_guided_min(struct subset * const ss);
int minimize_part_of_gm_set(struct subset *ss, int min_iterations);

struct subset *de_ss;
#endif /* __GUIDEDMIN_H__ */
