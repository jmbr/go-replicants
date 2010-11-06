#ifndef BITSET_H
#define BITSET_H
/**
 * @file bitset.h
 * @brief Interface to the bitset data structure.
 */

#include <stdlib.h>
#include <limits.h>


#define CHUNK_BITS	sizeof(int)*CHAR_BIT


/**
 * Bitset data structure.
 */
struct bitset {
	size_t length;          /**< Length in bits. */
	int chunk[];            /**< Storage for the set. */
};


extern struct bitset *new_bitset(size_t length);

extern void delete_bitset(struct bitset *self);

static inline void bitset_set(struct bitset *self, int index)
{
        div_t d = div(index, CHUNK_BITS);
        self->chunk[d.quot] |= (1 << d.rem);
}

static inline int bitset_get(const struct bitset *self, int index)
{
        div_t d = div(index, CHUNK_BITS);
        return (self->chunk[d.quot] >> d.rem) & 1;
}

static inline void bitset_clear(struct bitset *self, int index)
{
        div_t d = div(index, CHUNK_BITS);
        self->chunk[d.quot] &= ~(1 << d.rem);
}
#endif // !BITSET_H
