/**
 * @file test-protein.c
 * @brief Unit test for the protein data structure.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

#include "protein.h"
#include "utils.h"


static void test_constructor_and_destructor(void);
static void test_triple_scalar_product(void);

int main(void)
{
        test_constructor_and_destructor();
        test_triple_scalar_product();

        exit(EXIT_SUCCESS);
}



void test_constructor_and_destructor(void)
{
        const size_t num_atoms_1pgb = 56;  // Number of alpha carbons.
        const double positions_1pgb[][3] = {{ 13.935, 18.529, 29.843 }, { 13.088, 19.661, 26.283 }, { 12.726, 17.033, 23.612 }, { 12.179, 17.659, 19.887 }, { 10.253, 15.79, 17.221 }, { 11.082, 16.103, 13.475 }, { 8.009, 15.163, 11.389 }, { 8.628, 13.975, 7.913 }, { 5.213, 12.642, 6.966 }, { 3.589, 12.601, 3.497 }, { 1.291, 15.471, 4.538 }, { 2.624, 17.021, 7.787 }, { 6.24, 18.276, 7.903 }, { 8.146, 20.352, 10.511 }, { 9.429, 20.303, 14.17 }, { 7.799, 20.583, 17.589 }, { 9.057, 20.314, 21.095 }, { 7.744, 19.307, 24.468 }, { 8.907, 19.331, 28.095 }, { 8.985, 15.93, 29.75 }, { 10.637, 14.408, 32.787 }, { 12.217, 11.676, 30.658 }, { 12.452, 10.118, 27.226 }, { 9.683, 7.537, 27.584 }, { 7.249, 10.206, 28.617 }, { 8.146, 12.442, 25.683 }, { 7.82, 9.45, 23.25 }, { 4.309, 8.916, 24.682 }, { 3.067, 12.375, 23.977 }, { 4.843, 12.494, 20.533 }, { 3.247, 9.234, 19.455 }, { -0.277, 10.493, 20.449 }, { 0.341, 13.738, 18.651 }, { 1.53, 11.949, 15.49 }, { -1.401, 9.62, 15.646 }, { -3.917, 12.441, 16.16 }, { -2.542, 14.061, 13.009 }, { -2.684, 10.787, 11.031 }, { 1.109, 10.312, 10.729 }, { 2.165, 6.704, 10.736 }, { 5.826, 6.025, 10.072 }, { 9.219, 4.653, 11.307 }, { 10.629, 6.26, 14.43 }, { 14.117, 6.981, 15.677 }, { 15.402, 8.494, 18.951 }, { 18.798, 10.067, 19.525 }, { 19.598, 10.679, 23.19 }, { 22.584, 12.858, 22.345 }, { 20.338, 15.595, 20.936 }, { 17.089, 14.516, 22.58 }, { 15.427, 14.336, 19.205 }, { 12.767, 11.942, 17.92 }, { 11.971, 11.553, 14.286 }, { 8.973, 10.04, 12.475 }, { 9.147, 9.339, 8.722 }, { 6.283, 8.177, 6.48 }, };

        printf("%s\n", __FUNCTION__);

        struct protein *m = new_protein_1pgb();
        assert(m != NULL);

        assert(m->num_atoms == num_atoms_1pgb);

        print_matrix(stdout, m->positions);

        for (size_t i = 0; i < m->num_atoms; i++) {
                gsl_vector_view v = protein_get_atom(m, i);
                assert(gsl_vector_get(&v.vector, 0) == positions_1pgb[i][0]);
                assert(gsl_vector_get(&v.vector, 1) == positions_1pgb[i][1]);
                assert(gsl_vector_get(&v.vector, 2) == positions_1pgb[i][2]);
        }

        printf("%g\n", protein_distance(m, 0, 0));
        printf("%g\n", protein_distance(m, 0, 1));
        printf("%g\n", protein_distance(m, 1, 2));
        printf("%g\n", protein_distance(m, 0, 3));

        struct protein *g = new_protein_2gb1();
        assert(g != NULL);

        protein_plot(m);
        /* protein_plot(g); */

        delete_protein(m);
        delete_protein(g);
}



void test_triple_scalar_product(void)
{
        struct protein *m = new_protein_1pgb();
        assert(m != NULL);

        printf("%s\n", __FUNCTION__);

        gsl_vector_view u = protein_get_atom(m, 0);
        gsl_vector_view v = protein_get_atom(m, 1);
        gsl_vector_view w = protein_get_atom(m, 2);

        gsl_vector_fprintf(stdout, &u.vector, "u %g");
        gsl_vector_fprintf(stdout, &v.vector, "v %g");
        gsl_vector_fprintf(stdout, &w.vector, "w %g");

        double tsp = triple_scalar_product(&u.vector, &v.vector, &w.vector);
        printf("<u, v x w> = %g\n", tsp);
        assert(gsl_fcmp(tsp, -111.89, 1e-4) == 0);

        printf("%g\n", tsp/fabs(tsp));

        delete_protein(m);
}
