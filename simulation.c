/**
 * @file simulation.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <gsl/gsl_rng.h>

#include "protein.h"
#include "contact-map.h"
#include "simulation.h"


static const char gnuplot_command_line[] = GNUPLOT_EXECUTABLE
                                           " -noraise -persist";



struct simulation *new_simulation(enum simulation_protein p,
                                  double d_max, double a)
{
        if (d_max <= 0.0 || a <= 0.0)
                return NULL;
        
        struct simulation *s;

        if ((s = calloc(1, sizeof(struct simulation))) == NULL)
                return NULL;

        /* Initialize the molecule to be simulated. */
        s->protein = p == SIMULATION_1PGB ? new_protein_1pgb()
                                          : new_protein_2gb1();
        if (s->protein == NULL)
                delete_simulation(s);

        /* Initialize native contact map. */
        s->d_max = d_max;
        if ((s->native_map = new_contact_map(s->protein, s->d_max)) == NULL)
                delete_simulation(s);

        /* Initialize the pseudo-random number generator. */
        s->rng = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_env_setup();
        gsl_rng_set(s->rng, gsl_rng_default_seed);
        
        s->a = a;

        if ((s->gnuplot = popen(gnuplot_command_line, "w")) == NULL)
                delete_simulation(s);

        return s;
}


void delete_simulation(struct simulation *self)
{
        assert(self);

        if (self->gnuplot != NULL)
                fclose(self->gnuplot);
        if (self->protein != NULL)
                delete_contact_map(self->native_map);
        if (self->rng != NULL)
                gsl_rng_free(self->rng);
        if (self->protein != NULL)
                delete_protein(self->protein);
        free(self);
}
