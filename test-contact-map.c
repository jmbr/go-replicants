/**
 * @file test-contact-map.c
 */


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "protein.h"
#include "contact-map.h"


int main(void)
{
        struct protein *p = new_protein_2gb1();
        const double d_max = 7.5;
        struct contact_map *c = new_contact_map(p, d_max);

        (void) contact_map_plot(c);

        delete_contact_map(c);
        delete_protein(p);

        exit(EXIT_SUCCESS);
}
