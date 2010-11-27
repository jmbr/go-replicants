#undef NDEBUG
#include "molecular-simulator.h"


static bool plot_results = false;

static void test_initialization_and_finalization(void);
static void test_triple_scalar_product(void);
static void test_movements(void);


int main(int argc, char __attribute__((unused)) *argv[])
{
        int opt;

        while ((opt = getopt(argc, argv, "gh")) != -1) {
                switch (opt) {
                case 'g':
                        plot_results = true;
                        break;
                case 'h':
                        fprintf(stderr, "Usage: %s [-g]\n", argv[0]);
                        exit(EXIT_SUCCESS);
                }
        }

        test_initialization_and_finalization();
        test_triple_scalar_product();
        test_movements();

        exit(EXIT_SUCCESS);
}



void test_initialization_and_finalization(void)
{
        const size_t num_atoms_1pgb = 56;  // Number of alpha carbons.
        const double atoms_1pgb[][3] = {{ 13.935, 18.529, 29.843 }, { 13.088, 19.661, 26.283 }, { 12.726, 17.033, 23.612 }, { 12.179, 17.659, 19.887 }, { 10.253, 15.79, 17.221 }, { 11.082, 16.103, 13.475 }, { 8.009, 15.163, 11.389 }, { 8.628, 13.975, 7.913 }, { 5.213, 12.642, 6.966 }, { 3.589, 12.601, 3.497 }, { 1.291, 15.471, 4.538 }, { 2.624, 17.021, 7.787 }, { 6.24, 18.276, 7.903 }, { 8.146, 20.352, 10.511 }, { 9.429, 20.303, 14.17 }, { 7.799, 20.583, 17.589 }, { 9.057, 20.314, 21.095 }, { 7.744, 19.307, 24.468 }, { 8.907, 19.331, 28.095 }, { 8.985, 15.93, 29.75 }, { 10.637, 14.408, 32.787 }, { 12.217, 11.676, 30.658 }, { 12.452, 10.118, 27.226 }, { 9.683, 7.537, 27.584 }, { 7.249, 10.206, 28.617 }, { 8.146, 12.442, 25.683 }, { 7.82, 9.45, 23.25 }, { 4.309, 8.916, 24.682 }, { 3.067, 12.375, 23.977 }, { 4.843, 12.494, 20.533 }, { 3.247, 9.234, 19.455 }, { -0.277, 10.493, 20.449 }, { 0.341, 13.738, 18.651 }, { 1.53, 11.949, 15.49 }, { -1.401, 9.62, 15.646 }, { -3.917, 12.441, 16.16 }, { -2.542, 14.061, 13.009 }, { -2.684, 10.787, 11.031 }, { 1.109, 10.312, 10.729 }, { 2.165, 6.704, 10.736 }, { 5.826, 6.025, 10.072 }, { 9.219, 4.653, 11.307 }, { 10.629, 6.26, 14.43 }, { 14.117, 6.981, 15.677 }, { 15.402, 8.494, 18.951 }, { 18.798, 10.067, 19.525 }, { 19.598, 10.679, 23.19 }, { 22.584, 12.858, 22.345 }, { 20.338, 15.595, 20.936 }, { 17.089, 14.516, 22.58 }, { 15.427, 14.336, 19.205 }, { 12.767, 11.942, 17.92 }, { 11.971, 11.553, 14.286 }, { 8.973, 10.04, 12.475 }, { 9.147, 9.339, 8.722 }, { 6.283, 8.177, 6.48 }, };

        assert(new_protein(0, (void *) 1) == NULL);
        assert(new_protein(1, NULL) == NULL);

        struct protein *m = new_protein_1pgb();
        assert(m != NULL);

        assert(m->num_atoms == num_atoms_1pgb);

        /* print_matrix(stdout, m->atoms); */

        for (size_t i = 0; i < m->num_atoms; i++) {
                assert(gsl_vector_get(m->atom[i], 0) == atoms_1pgb[i][0]);
                assert(gsl_vector_get(m->atom[i], 1) == atoms_1pgb[i][1]);
                assert(gsl_vector_get(m->atom[i], 2) == atoms_1pgb[i][2]);
        }

        printf("%g\n", protein_distance(m, 0, 0));
        printf("%g\n", protein_distance(m, 0, 1));
        printf("%g\n", protein_distance(m, 1, 2));
        printf("%g\n", protein_distance(m, 0, 3));

        struct protein *g = new_protein_2gb1();
        assert(g != NULL);

        if (plot_results) {
                FILE *gnuplot;

                gnuplot = popen("gnuplot -persist", "w");
                assert(gnuplot != NULL);
                protein_plot(m, gnuplot, "1PGB");
                pclose(gnuplot);

                gnuplot = popen("gnuplot -persist", "w");
                assert(gnuplot != NULL);
                protein_plot(g, gnuplot, "2GB1");
                pclose(gnuplot);
        }

        delete_protein(m);
        delete_protein(g);
}



void test_triple_scalar_product(void)
{
        struct protein *m = new_protein_1pgb();
        assert(m != NULL);

        gsl_vector *u = m->atom[0], *v = m->atom[1], *w = m->atom[2];

        print_vector(stdout, u);
        print_vector(stdout, v);
        print_vector(stdout, w);

        double tsp = triple_scalar_product(u, v, w);
        printf("<u, v x w> = %g\n", tsp);
        assert(gsl_fcmp(tsp, -111.89, 1e-4) == 0);

        printf("signum == %g\n", tsp/fabs(tsp));

        delete_protein(m);
}



void test_movements(void)
{
        if (!plot_results)
                return;

        /* const double atoms[][3] = {{-3.8, 0, 0}, */
        /*                            {0, 0, 1}, */
        /*                            {3.8, 0, 0}, */
        /*                            {2.0*3.8, 0, 0}, */
        /*                            {3.0*3.8, 0, 0}, */
        /*                            {4.0*3.8, 0, 0}, */
        /*                            {5.0*3.8, 0, 0}}; */
        /* const size_t num_atoms = 7; */

        bool undo = true;
        struct protein *n;
        struct protein *(*constructor)(void) = new_protein_1pgb;

        FILE *g = popen("gnuplot -persist", "w");
        assert(g != NULL);

        gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_env_setup();
        gsl_rng_set(r, (unsigned) time(NULL));

        n = constructor();
        for (size_t i = 0; i < 100; i++) {
                protein_plot(n, g, "End chain movement (1)");
                protein_do_movement(n, r, PROTEIN_END_MOVE_FIRST, 0, undo);
        }
        delete_protein(n);
        sleep(3);

        n = constructor();
        for (size_t i = 0; i < 100; i++) {
                protein_do_movement(n, r, PROTEIN_END_MOVE_LAST, 0, undo);
                protein_plot(n, g, "End chain movement (2)");
        }
        delete_protein(n);
        sleep(3);

        n = constructor();
        for (size_t i = 0; i < 250; i++) {
                protein_do_movement(n, r, PROTEIN_SHIFT_MOVE, 1, undo);
                protein_plot(n, g, "Shifting movement");
        }
        delete_protein(n);
        sleep(3);

        n = constructor();
        for (size_t i = 0; i < 250; i++) {
                protein_do_movement(n, r, PROTEIN_SPIKE_MOVE, 1, undo);
                protein_plot(n, g, "Spike movement");
        }
        delete_protein(n);
        sleep(3);

        n = constructor();
        for (size_t i = 0; i < 250; i++) {
                protein_do_movement(n, r, PROTEIN_PIVOT_MOVE, 52, undo);
                protein_plot(n, g, "Pivoting movement");
        }
        delete_protein(n);

        n = constructor();
        protein_scramble(n, r);
        protein_plot(n, g, "Scrambled molecule");
        delete_protein(n);

        pclose(g);

        gsl_rng_free(r);
}
