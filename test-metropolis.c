/**
 * @file test-metropolis.c
 * @brief Unit test for the implementation of the Metropolis-Hastings algorithm.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

#include "metropolis.h"


const double ratio_of_zeros_to_ones = 0.1;


static metropolis_state
generate_state(gsl_rng *rng,
               metropolis_state __attribute__((unused)) previous)
{
        return (metropolis_state) ((int) gsl_rng_uniform_int(rng, 2));
}

static double
pdf_factor(metropolis_state x)
{
        return ((int) x == 0) ? ratio_of_zeros_to_ones
                              : (1 - ratio_of_zeros_to_ones);
}

int main(void)
{
        unsigned long output[2] = {0, 0};
        const int num_iters = 1000000;
        metropolis_state s = (metropolis_state) 0;
        gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, (unsigned) time(NULL));

        struct metropolis *m = new_metropolis(rng, generate_state, pdf_factor);

        for (int i = 0; i < num_iters; i++) {
                s = metropolis_iteration(m, s);
                ++output[(int) s];
        }

        printf("# of zeros = %lu, # of ones = %lu\n", output[0], output[1]);
        printf("ratio of zeros to ones = %g\n", (double) output[0]/num_iters);
        printf("Acceptance ratio = %g\n", metropolis_get_acceptance_ratio(m));

        gsl_rng_free(rng);
        delete_metropolis(m);
        return gsl_fcmp((double) output[0]/(double) num_iters,
                        ratio_of_zeros_to_ones, 1e-2);
}
