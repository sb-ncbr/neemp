/*
 * NEEMP - bitarray.c
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "bitarray.h"
#include "neemp.h"

static inline int b_idx(int bitno);
static inline int b_off(int bitno);

/* Initialize bit array */
void b_init(struct bit_array * const b, int count) {

	assert(b != NULL);

	b->ints_count = count / sizeof(int) + 1;
	b->bits = (int *) calloc(b->ints_count, sizeof(int));
	if(!b->bits)
		EXIT_ERROR(MEM_ERROR, "%s", "Cannot allocate memory for bitarray.\n");
}

/* Set bits in bdest to be the same as in bsrc */
void b_set_as(struct bit_array * const bdest, const struct bit_array * const bsrc) {

	assert(bdest != NULL);
	assert(bsrc != NULL);
	assert(bdest->bits_count == bsrc->bits_count);
	assert(bdest->ints_count == bsrc->ints_count);

	memcpy(bdest->bits, bsrc->bits, bsrc->ints_count * sizeof(int));
}

/* Destroy contents of bit array */
void b_destroy(struct bit_array * const b) {

	assert(b != NULL);

	free(b->bits);
}

/* Return number of ints befor bitno */
static inline int b_idx(int bitno) {

	return bitno / sizeof(int);
}

/* Return offset to the int containing bitno */
static inline int b_off(int bitno) {

	return bitno % sizeof(int);
}

/* Return the value of bit bitno */
int b_get(const struct bit_array * const b, int bitno) {

	assert(b != NULL);
	assert(bitno < b->bits_count);

	return b->bits[b_idx(bitno)] & (1 << b_off(bitno)) ? 1 : 0;
}

/* Set bit bitno */
void b_set(struct bit_array * const b, int bitno) {

	assert(b != NULL);
	assert(bitno < b->bits_count);

	b->bits[b_idx(bitno)] |= 1 << b_off(bitno);
}

/* Clear bit bitno */
void b_clear(struct bit_array * const b, int bitno) {

	assert(b != NULL);
	assert(bitno < b->bits_count);

	b->bits[b_idx(bitno)] &= ~(1 << b_off(bitno));
}

/* Flip bit bitno, return previous value */
int b_flip(struct bit_array * const b, int bitno) {

	assert(b != NULL);
	assert(bitno < b->bits_count);

	int previous = b_get(b, bitno);

	b->bits[b_idx(bitno)] ^= 1 << b_off(bitno);

	return previous;
}

/* Set all bit to 1 */
void b_set_all(struct bit_array * const b) {

	assert(b != NULL);

	memset(b->bits, 0xff, b->ints_count * sizeof(int));
}

/* Set all bit to 0 */
void b_clear_all(struct bit_array * const b) {

	assert(b != NULL);

	memset(b->bits, 0x0, b->ints_count * sizeof(int));
}

/* Return number of bits set */
int b_count_bits(const struct bit_array * const b) {

	assert(b != NULL);

	int count = 0;
	for(int i = 0; i < b->bits_count; i++)
		if(b_get(b, i))
			count++;

	return count;
}
