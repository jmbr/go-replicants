/*
 * test-metropolis.c -- Unit tests for the implementation of the Metropolis-Hastings algorithm.
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <math.h>
#include <time.h>

#include "metropolis.h"


static state gen(const struct metropolis *m, state previous)
{
        return (state) (((int) round(metropolis_random(m)*10)) % 2);
}

static double pi(state x)
{
        return ((int) x == 0) ? 5.0 : 5.0;
}

int main(int argc, char *argv[])
{
        struct metropolis *m = new_metropolis(time(NULL), gen, pi);
        int output[2] = {0, 0};
        const unsigned num_iters = 1e6;
        state s = (state) 0;
        for (int i = 0; i < num_iters; i++) {
                s = metropolis_iteration(m, s);
                ++output[(int) s];
        }

        double ratio = (double) output[0]/num_iters;
        printf("Ratio: %g ... %s.\n", ratio, (fabs(ratio - 0.5) < 1e-3) ? "PASS" : "FAIL");

        delete_metropolis(m);
        exit(EXIT_SUCCESS);
}
