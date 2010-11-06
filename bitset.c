/**
 * @file bitset.c
 * @brief Implementation of a bit set data structure.
 */

#include <stdlib.h>
#include <assert.h>

#include "bitset.h"


static inline unsigned int num_chunks(size_t length)
{
        div_t d = div((int) length, CHUNK_BITS);
        return (unsigned int) (d.rem == 0 ? d.quot : d.quot+1);
}

struct bitset *new_bitset(size_t length)
{
	struct bitset *b;
        size_t total_size;

        if (length == 0)
                return NULL;

        total_size = sizeof(struct bitset) + num_chunks(length)*sizeof(int);

        if ((b = calloc(1, total_size)) == NULL)
		return NULL;

	b->length = length;

	return b;
}

void delete_bitset(struct bitset *self)
{
	assert(self);

	free(self);
}
