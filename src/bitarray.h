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

#ifndef __BITARRAY_H__
#define __BITARRAY_H__

struct bit_array {

	int *bits;
	int bits_count;
	int ints_count;
};

void b_init(struct bit_array * const b, int count);
void b_set_as(struct bit_array * const bdest, const struct bit_array * const bsrc);
void b_destroy(struct bit_array * const b);
int b_get(const struct bit_array * const b, int bitno);
void b_set(struct bit_array * const b, int bitno);
void b_clear(struct bit_array * const b, int bitno);
int b_flip(struct bit_array * const b, int bitno);
void b_set_all(struct bit_array * const b);
void b_clear_all(struct bit_array * const b);
int b_count_bits(const struct bit_array * const b);

#endif /* __BITARRAY_H__ */
