#undef NDEBUG
#include "molecular-simulator.h"


static void show_progress(struct replicas *r, size_t k);


int main(void)
{

        struct protein *p = new_protein_1pgb();
        assert(p != NULL);

        /* Initialize the pseudo-random number generator. */
        gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
        assert(rng != NULL);
        gsl_rng_env_setup();
        gsl_rng_set(rng, gsl_rng_default_seed);

        const double alpha0 = 0.0001;
        const size_t num_replicas = 8;
        double temperatures[num_replicas];
        linspace(0.05, 0.9, num_replicas, temperatures);
        struct simulation_options options = { .d_max = 10.0, .a = 0.5 };
        struct replicas *r = new_replicas(p, rng, alpha0,
                                          num_replicas, temperatures,
                                          &options);

        size_t k;
        for (replicas_first_iteration(r), k = 1;
             replicas_have_not_converged(r);
             replicas_next_iteration(r), k++)
        {
                show_progress(r, k);
        }

        delete_protein(p);
        delete_replicas(r);
        gsl_rng_free(rng);
        exit(EXIT_SUCCESS);
}



void show_progress(struct replicas *r, size_t k)
{
        if (k != 1 && k % 10000 != 0)
                return;

        printf("total number of exchanged replicas: %u.\n",
               replicas_total_exchanges(r));
        for (size_t j = 0; j < r->num_replicas; j++) {
                struct simulation *s = r->replica[j];
                simulation_print_info(s, stdout);
                protein_plot(s->protein, s->gnuplot,
                             "Replica %u, iteration %u (%g).  Acceptance: %2.03g", j,
                             k, s->energy, simulation_get_acceptance_ratio(s));
                printf("progress for replica %u: %1.3g%% (%g/%g)\n", j,
                       100*s->energy/s->orig_energy, s->energy, s->orig_energy);
        }
        fflush(stdout);
}
