/*
 * metropolis.c -- Implementation of the Metropolis-Hastings algorithm.
 */

#include <assert.h>
#include <gsl/gsl_rng.h>

#include "metropolis.h"


struct metropolis {
        gsl_rng *rng;
        generator_fn generator;
        pi_fn pi;
};

struct metropolis *new_metropolis(unsigned int seed,
                                  generator_fn generator,
                                  pi_fn pi)
{
        struct metropolis *m;

        m = malloc(sizeof(struct metropolis)); /* XXX */
        m->rng = gsl_rng_alloc(gsl_rng_mt19937);
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

state metropolis_iteration(const struct metropolis *self, state current_state)
{
        assert(self);

        state candidate = self->generator(self, current_state);

        return (metropolis_random(self) < self->pi(candidate)/self->pi(current_state))
                ? candidate
                : current_state;
}
