#undef NDEBUG
#include "molecular-simulator.h"


static void show_progress(struct replicas *r, size_t k);


int main(void)
{

        struct protein *p = new_protein_2gb1();
        assert(p != NULL);

        /* Initialize the pseudo-random number generator. */
        gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
        assert(rng != NULL);
        gsl_rng_env_setup();
        gsl_rng_set(rng, gsl_rng_default_seed);

        double temperatures[] = {0.100, 0.105, 0.111, 0.118, 0.126, 0.134, 0.144, 0.156, 0.169, 0.185, 0.204, 0.228, 0.258, 0.298, 0.351, 0.428, 0.547, 0.760, 1.242, };
        const size_t num_replicas = sizeof(temperatures)/sizeof(temperatures[0]);
        /* linspace(0.1, 0.9, num_replicas, temperatures); */
        struct simulation_options options = {
                .rng = rng, .a = 0.5, .d_max = 10.0,
                .num_replicas = num_replicas, .temperatures = temperatures
        };
        struct replicas *r = new_replicas(p, &options);

        size_t k = 1;
        replicas_thermalize(r, 1);
        while (replicas_have_not_converged(r)) {
                replicas_next_iteration(r);
                show_progress(r, k);
                ++k;
        }

        delete_protein(p);
        delete_replicas(r);
        gsl_rng_free(rng);
        exit(EXIT_SUCCESS);
}



void show_progress(struct replicas *r, size_t k)
{
        if (k != 1 && k % 100 != 0)
                return;

        printf("total number of exchanges: %u\n", replicas_total_exchanges(r));
        
        replicas_print_info(r, stdout);

        for (size_t j = 0; j < r->num_replicas; j++) {
                struct simulation *s = r->replica[j];

                simulation_print_info(s, stdout);

                /* protein_plot(s->protein, s->gnuplot, */
                /*              "Replica %u, iteration %u (%g).  Acceptance: %2.03g", j, */
                /*              k, s->energy, simulation_get_acceptance_ratio(s)); */

                printf("progress for replica %u: %g/%g\n",
                       j, s->energy, r->orig_energy);
        }
        fflush(stdout);
}
