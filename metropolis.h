#ifndef METROPOLIS_H
#define METROPOLIS_H
/**
 * @file metropolis.h
 * @brief Interface to the Metropolis-Hastings algorithm.
 */

#include <stdint.h>

#include <gsl/gsl_rng.h>


typedef void *metropolis_state;
typedef metropolis_state (*metropolis_generator)(gsl_rng *, metropolis_state);
typedef double (*metropolis_pdf_factor)(metropolis_state);
/**
 * Metropolis structure.
 * Provides context for the Metropolis-Hastings algorithm.
 */
struct metropolis {
        /** Random number generator. */
        gsl_rng *rng;

        /** Pointer to a function that generates new states. */
        metropolis_generator generator;

        /** Pointer to a function proportional to the p.d.f. of the
         * target probability distribution. */
        metropolis_pdf_factor pdf_factor;

        /** Number of accepted moves and number of total moves. */
        uint_fast32_t accepted, total;
};


extern struct metropolis *new_metropolis(gsl_rng *rng,
                                         metropolis_generator generator,
                                         metropolis_pdf_factor pdf_factor);

extern void delete_metropolis(struct metropolis *self);

extern metropolis_state metropolis_iteration(struct metropolis *self,
                                             metropolis_state current_state);

extern double metropolis_get_acceptance_ratio(const struct metropolis *self);
#endif // !METROPOLIS_H
