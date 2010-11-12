/**
 * @file test-contact-map.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <gsl/gsl_math.h>

#include "protein.h"
#include "potential.h"
#include "contact-map.h"


static void test_initialization_and_finalization(void);
static void test_potential_energy(void);


int main(void)
{
        test_initialization_and_finalization();
        test_potential_energy();

        exit(EXIT_SUCCESS);
}


void test_initialization_and_finalization(void)
{
        struct protein *p = new_protein_2gb1();
        assert(p != NULL);

        struct contact_map *c1 = new_contact_map(p, 10);
        struct contact_map *c2 = new_contact_map(p, 7.5);
        assert(c1 != NULL && c2 != NULL);
        
        printf("The following plots should be compared against those appearing"
               "in Figure 3.2 of Lidia Prieto's Ph.D. thesis (p. 60).\n");
        
        contact_map_plot(c1);
        contact_map_plot(c2);

        delete_contact_map(c1);
        delete_contact_map(c2);
        delete_protein(p);
}

void test_potential_energy(void)
{
        struct protein *p1 = new_protein_1pgb();
        struct protein *p2 = new_protein_2gb1();
        assert(p1 != NULL && p2 != NULL);
        
        const double d_max = 10;
        struct contact_map *c1 = new_contact_map(p1, d_max);
        struct contact_map *c2 = new_contact_map(p2, d_max);
        assert(c1 != NULL && c2 != NULL);

        double p_1pgb = potential(p1, c1, 0.9);
        double p_2gb1 = potential(p2, c2, 0.9);

        printf("U(1pgb) = %g\n", p_1pgb);
        printf("U(2gb1) = %g\n", p_2gb1);

        assert(gsl_fcmp(p_1pgb, -460.0, 1e-15) == 0);
        assert(gsl_fcmp(p_2gb1, -460.0, 1e-15) == 0);

        delete_contact_map(c1);
        delete_contact_map(c2);
        delete_protein(p1);
        delete_protein(p2);
}
