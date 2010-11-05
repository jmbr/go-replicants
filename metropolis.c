/**
 * @file metropolis.c
 * @brief Implementation of the Metropolis-Hastings algorithm.
 */

#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_rng.h>

#include "metropolis.h"


struct metropolis {
        gsl_rng *rng;
        metropolis_generator generator;
        metropolis_pi pi;
};


struct metropolis *new_metropolis(unsigned int seed,
                                  metropolis_generator generator,
                                  metropolis_pi pi)
{
        struct metropolis *m;

        if ((m = malloc(sizeof(struct metropolis))) == NULL)
                return NULL;

        if ((m->rng = gsl_rng_alloc(gsl_rng_mt19937)) == NULL) {
                free(m);
                return NULL;
        }

        gsl_rng_set(m->rng, seed);
        m->generator = generator;
        m->pi = pi;

        return m;
}

void delete_metropolis(struct metropolis *self)
{
        assert(self);
        assert(self->rng);

        gsl_rng_free(self->rng);
        free(self);
}

double metropolis_random(const struct metropolis *self)
{
        assert(self);
        assert(self->rng);

        return gsl_rng_uniform(self->rng);
}

metropolis_state metropolis_iteration(const struct metropolis *self, metropolis_state current_state)
{
        assert(self);

        metropolis_state candidate = self->generator(current_state);

        return (metropolis_random(self) < self->pi(candidate)/self->pi(current_state))
                ? candidate
                : current_state;
}
