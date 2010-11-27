#include "molecular-simulator.h"


struct replicas *new_replicas(const struct protein *p,
                              gsl_rng *rng, double alpha0,
                              size_t num_replicas, double temperatures[],
                              const struct simulation_options *options)
{
        if (p == NULL || rng == NULL
            || alpha0 < 0.0 || alpha0 > 1.0
            || num_replicas == 0 || options == NULL)
        {
                return NULL;
        }
        
        struct replicas *r;
        r = calloc(1, sizeof(struct replicas)
                      + num_replicas*sizeof(struct simulation *));
        if (r == NULL)
                return NULL;

        r->rng = rng;

        r->alpha0 = alpha0;

        r->exchanges = calloc(num_replicas, sizeof(uint_fast32_t));
        r->total = 0;

        r->num_replicas = num_replicas;

        for (size_t k = 0; k < r->num_replicas; k++)
                r->replica[k] = new_simulation(p, rng,
                                               temperatures[k], options);

        return r;
}

void delete_replicas(struct replicas *self)
{
        for (size_t k = 0; k < self->num_replicas; k++)
                delete_simulation(self->replica[k]);
        free(self->exchanges);
        free(self);
}



void replicas_first_iteration(struct replicas *self)
{
        size_t k;
        #pragma omp parallel for private(k)
        for (k = 0; k < self->num_replicas; k++)
                simulation_first_iteration(self->replica[k]);
}

void replicas_next_iteration(struct replicas *self)
{
        double u = gsl_rng_uniform(self->rng);

        dprintf("u == %g, alpha0 == %g\n", u, self->alpha0);

        if (u > self->alpha0) {
                size_t k;
                #pragma omp parallel for private(k)
                for (k = 0; k < self->num_replicas; k++)
                        simulation_next_iteration(self->replica[k]);
        } else {
                ++self->total;

                size_t k = gsl_rng_uniform_int(self->rng,
                                               self->num_replicas-1);

                struct simulation *s1 = self->replica[k];
                struct simulation *s2 = self->replica[k+1];
                const double U1 = potential(s1->protein,
                                            s1->native_map, s1->a);
                const double U2 = potential(s2->protein,
                                            s2->native_map, s2->a);
                const double DU = U2 - U1;
                const double DB = 1.0/s2->T - 1.0/s1->T;
                const double p = exp(DB*DU);

                const double r = gsl_rng_uniform(self->rng);
                if (r < gsl_min(1.0, p)) {
                        dprintf("swapping replicas %u and %u.\n", k, k+1);
                        struct protein *x = s1->protein;
                        s1->protein = s2->protein;
                        s1->energy = U2;
                        s2->protein = x;
                        s2->energy = U1;
                        ++self->exchanges[k];
                }
        }
}



bool replicas_have_converged(const struct replicas *self)
{
        bool status = false;

        for (size_t s = 0; s < self->num_replicas; s++)
                status |= simulation_has_converged(self->replica[s]);

        return status;
}



uint_fast32_t replicas_total_exchanges(const struct replicas *self)
{
        uint_fast32_t total = 0;
        for (size_t r = 0; r < self->num_replicas - 1; r++)
                total += self->exchanges[r];
        return total;
}
