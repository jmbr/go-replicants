/*
 * test-metropolis.c -- Unit tests for the implementation of the Metropolis-Hastings algorithm.
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>

#include "metropolis.h"


static gsl_rng *rng;


static state gen(state previous)
{
        return (state) (((int) round(gsl_rng_uniform(rng)*10)) % 2);
}

static double pi(state x)
{
        return ((int) x == 0) ? 1.0 : 9.0;
}

int main(int argc, char *argv[])
{
        struct metropolis *m = new_metropolis(time(NULL), gen, pi);
        double output[2] = {0, 0};
        const unsigned num_iters = 1e7;
        state s = (state) 0;

        rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, time(NULL));

        for (int i = 0; i < num_iters; i++) {
                s = metropolis_iteration(m, s);
                ++output[(int) s];
        }

        delete_metropolis(m);
        exit((fabs(output[0]/num_iters - 0.1) < 1e-2) ? EXIT_SUCCESS : EXIT_FAILURE);
}
