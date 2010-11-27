#undef NDEBUG
#include "molecular-simulator.h"


static bool plot_results = false;


static struct protein *new_truncated_1pgb(size_t start, size_t end);

static void test_simulation(void);


int main(int argc, char *argv[])
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

        test_simulation();

        exit(EXIT_SUCCESS);
}



static void show_progress(struct simulation *s, size_t k)
{
        if (k != 1 || k % 10000 != 0)
                return;

        if (plot_results)
                protein_plot(s->protein, s->gnuplot,
                             "#%u (%g).  Acceptance: %2.03g%%",
                             k, s->energy, simulation_get_acceptance_ratio(s));

        printf("Progress %1.3g%% (%g/%g)\n",
               100*s->energy/s->orig_energy, s->energy, s->orig_energy);
        fflush(stdout);
}

static void simulate_fragment(size_t start, size_t end)
{
        gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_env_setup();
        gsl_rng_set(rng, gsl_rng_default_seed);

        struct simulation_options opts = { .d_max = 10, .a = 0.7 };
        struct protein *p = new_truncated_1pgb(start, end);
        struct simulation *s = new_simulation(p, rng, 1.0, &opts);
        assert(p != NULL && s != NULL);
        delete_protein(p);

        printf("Simulating fragment #%2u-%2u of 1PGB.\n", start, end);

        printf("fragment = [\n");
        for (size_t k = 0; k < s->protein->num_atoms; k++)
                print_vector(stdout, s->protein->atom[k]);
        printf("]\n");

        if (plot_results) {
                protein_plot(s->protein, s->gnuplot, "Original: (%3.3g)",
                             s->energy);
                /* contact_map_plot(s->native_map, s->gnuplot); */
                sleep(3);
        }

        size_t k = 0;
        for (simulation_first_iteration(s);
             simulation_has_not_converged(s);
             simulation_next_iteration(s), k++)
        {
                show_progress(s, k);
        }

        printf("\nEnergy after simulation: %g\n", s->energy);

        gsl_rng_free(rng);
        delete_simulation(s);
}

void test_simulation(void)
{
        simulate_fragment(0, 3);
        simulate_fragment(6, 10);
        /* simulate_fragment(0, 15); */
        /* simulate_fragment(20, 30); */
}

struct protein *new_truncated_1pgb(size_t start, size_t end)
{
        const double positions[][3] = {
                { 13.935, 18.529, 29.843 },
                { 13.088, 19.661, 26.283 },
                { 12.726, 17.033, 23.612 },
                { 12.179, 17.659, 19.887 },
                { 10.253, 15.79, 17.221 },
                { 11.082, 16.103, 13.475 },
                { 8.009, 15.163, 11.389 },
                { 8.628, 13.975, 7.913 },
                { 5.213, 12.642, 6.966 },
                { 3.589, 12.601, 3.497 },
                { 1.291, 15.471, 4.538 },
                { 2.624, 17.021, 7.787 },
                { 6.24, 18.276, 7.903 },
                { 8.146, 20.352, 10.511 },
                { 9.429, 20.303, 14.17 },
                { 7.799, 20.583, 17.589 },
                { 9.057, 20.314, 21.095 },
                { 7.744, 19.307, 24.468 },
                { 8.907, 19.331, 28.095 },
                { 8.985, 15.93, 29.75 },
                { 10.637, 14.408, 32.787 },

                { 12.217, 11.676, 30.658 },
                { 12.452, 10.118, 27.226 },
                { 9.683, 7.537, 27.584 },
                { 7.249, 10.206, 28.617 },
                { 8.146, 12.442, 25.683 },
                { 7.82, 9.45, 23.25 },
                { 4.309, 8.916, 24.682 },
                { 3.067, 12.375, 23.977 },
                { 4.843, 12.494, 20.533 },
                { 3.247, 9.234, 19.455 },
                { -0.277, 10.493, 20.449 },
                { 0.341, 13.738, 18.651 },
                { 1.53, 11.949, 15.49 },
                { -1.401, 9.62, 15.646 },
                { -3.917, 12.441, 16.16 },
                { -2.542, 14.061, 13.009 },
                { -2.684, 10.787, 11.031 },

                { 1.109, 10.312, 10.729 },
                { 2.165, 6.704, 10.736 },
                { 5.826, 6.025, 10.072 },
                { 9.219, 4.653, 11.307 },
                { 10.629, 6.26, 14.43 },
                { 14.117, 6.981, 15.677 },
                { 15.402, 8.494, 18.951 },
                { 18.798, 10.067, 19.525 },
                { 19.598, 10.679, 23.19 },
                { 22.584, 12.858, 22.345 },
                { 20.338, 15.595, 20.936 },
                { 17.089, 14.516, 22.58 },
                { 15.427, 14.336, 19.205 },
                { 12.767, 11.942, 17.92 },
                { 11.971, 11.553, 14.286 },
                { 8.973, 10.04, 12.475 },
                { 9.147, 9.339, 8.722 },
                { 6.283, 8.177, 6.48 },
        };
        const size_t num_atoms = sizeof(positions)/sizeof(positions[0]);

        if (start > end || end >= num_atoms)
                return NULL;
        else
                return new_protein(end - start + 1,
                                   (const double *) &positions[start]);
}
