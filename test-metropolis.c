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


static gsl_rng *rng;

const double ratio_of_zeros_to_ones = 0.1;


static metropolis_state gen(metropolis_state __attribute__((unused)) previous)
{
        return (metropolis_state) (((int) round(gsl_rng_uniform(rng)*10)) % 2);
}

static double pi(metropolis_state x)
{
        return ((int) x == 0) ? ratio_of_zeros_to_ones : (1 - ratio_of_zeros_to_ones);
}

int main(void)
{
        unsigned long output[2] = {0, 0};
        const int num_iters = 1000000;
        metropolis_state s = (metropolis_state) 0;
        unsigned seed = (unsigned) time(NULL);
        struct metropolis *m = new_metropolis(seed, gen, pi);

        rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, seed);

        for (int i = 0; i < num_iters; i++) {
                s = metropolis_iteration(m, s);
                ++output[(int) s];
        }

        printf("# of zeros = %lu, # of ones = %lu\n", output[0], output[1]);
        printf("Acceptance ratio = %g\n", metropolis_get_acceptance_ratio(m));

        gsl_rng_free(rng);
        delete_metropolis(m);
        return gsl_fcmp((double) output[0]/num_iters, ratio_of_zeros_to_ones, 0.1);
}
