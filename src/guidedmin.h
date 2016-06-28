/* Copyright 2013-2016 Jana Pazurikova (pazurikova@ics.muni.cz)
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

#ifndef __GUIDEDMIN_H__
#define __GUIDEDMIN_H__

#include "subset.h"
#include "diffevolution.h"

void run_guided_min(struct subset * const ss);
int minimize_part_of_gm_set(struct subset *ss, int min_iterations);

struct subset *de_ss;
#endif /* __GUIDEDMIN_H__ */
