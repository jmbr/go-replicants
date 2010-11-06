/**
 * @file test-bitset.c
 * @brief Unit test for the bitset data structure.
 */

#include <stdio.h>
#include <assert.h>

#include "bitset.h"


int main(void)
{
        struct bitset *b = new_bitset(1);
        assert(b != NULL);

        assert(bitset_get(b, 0) == 0);
        bitset_set(b, 0);
        assert(bitset_get(b, 0) == 1);

        delete_bitset(b);



        const int N = 100;
        struct bitset *a = new_bitset((size_t) N);

        for (int i = 0; i < N; i++) {
                if (i % 3 != 0) {
                        bitset_set(a, i);
                        putchar('1');
                } else {
                        putchar('0');
                }
        }
        putchar('\n');
        
        for (int i = 0; i < N; i++) {
                int s = bitset_get(a, i);
                if (i % 3 == 0)
                        assert(s == 0);
                else
                        assert(s == 1);
                putchar(s == 0 ? '0' : '1');
        }
        putchar('\n');
        
        delete_bitset(a);

        return 0;
}
