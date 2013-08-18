/*
 * NEEMP - bitarray.h
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

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
