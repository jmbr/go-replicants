/**
 * @file metropolis.c
 * @brief Implementation of the Metropolis-Hastings algorithm.
 */

#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_rng.h>

#include "metropolis.h"


struct metropolis *new_metropolis(gsl_rng *rng,
                                  metropolis_generator generator,
                                  metropolis_pdf_factor pdf_factor)
{
        if (rng == NULL)
                return NULL;

        struct metropolis *m = calloc(1, sizeof(struct metropolis));
        if (m == NULL)
                return NULL;

        m->rng = rng;
        m->generator = generator;
        m->pdf_factor = pdf_factor;
        m->accepted = 0;
        m->total = 0;

        return m;
}

void delete_metropolis(struct metropolis *self)
{
        assert(self);

        free(self);
}

metropolis_state metropolis_iteration(struct metropolis *self,
                                      metropolis_state current_state)
{
        assert(self != NULL);

        ++self->total;

        metropolis_state candidate_state = self->generator(self->rng,
                                                           current_state);
        const double p_current = self->pdf_factor(current_state);
        const double p_candidate = self->pdf_factor(candidate_state);

        if (gsl_rng_uniform(self->rng) < p_candidate/p_current) {
                ++self->accepted;
                return candidate_state;
        } else
                return current_state;
}

double metropolis_get_acceptance_ratio(const struct metropolis *self)
{
        assert(self != NULL);

        return (double) self->accepted/(double) self->total;
}
