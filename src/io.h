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

#ifndef __IO_H__
#define __IO_H__

#include "subset.h"

void load_molecules(void);
void load_charges(void);
void load_parameters(struct kappa_data * const ss);
void load_user_atom_types(void);

void output_charges(const struct subset * const ss);
void output_charges_stats(const struct subset * const ss);
void output_parameters(const struct subset * const ss);

#endif /* __IO_H__ */
