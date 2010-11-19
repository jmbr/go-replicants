/**
 * @file simulation.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_const.h>
#include <gsl/gsl_rng.h>

#include "protein.h"
#include "potential.h"
#include "contact-map.h"
#include "simulation.h"


const char gnuplot_command_line[] = GNUPLOT_EXECUTABLE " -persist";



struct simulation *new_simulation(struct protein *p,
                                  double d_max, double a)
{
        if (p == NULL || d_max <= 0.0 || a <= 0.0)
                return NULL;

        struct simulation *s;

        if ((s = calloc(1, sizeof(struct simulation))) == NULL)
                return NULL;

        s->protein = p;

        /* Initialize native contact map. */
        s->d_max = d_max;
        if ((s->native_map = new_contact_map(s->protein, s->d_max)) == NULL)
                delete_simulation(s);

        s->a = a;

        /* Compute the potential energy for the protein. */
        s->U = potential(s->protein, s->native_map, s->a);

        /* Initialize the pseudo-random number generator. */
        s->rng = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_env_setup();
        gsl_rng_set(s->rng, gsl_rng_default_seed);

        if ((s->gnuplot = popen(gnuplot_command_line, "w")) == NULL)
                delete_simulation(s);

        s->accepted = 0;
        s->total = 0;

        s->configurations = fopen("configurations.dat", "w");
        s->energies = fopen("energies.dat", "w");
        /* XXX Error-checking. */

        return s;
}

void delete_simulation(struct simulation *self)
{
        assert(self);

        if (self->gnuplot != NULL)
                pclose(self->gnuplot);
        if (self->protein != NULL)
                delete_contact_map(self->native_map);
        if (self->rng != NULL)
                gsl_rng_free(self->rng);
        if (self->configurations != NULL)
                fclose(self->configurations);
        if (self->energies != NULL)
                fclose(self->energies);

        free(self);
}



double simulation_get_acceptance_ratio(const struct simulation *self)
{
        assert(self != NULL);

        return (double) self->accepted/(double) self->total;
}



static inline double
compute_potential_energy(const struct protein *x, const struct simulation *s)
{
        return potential(x, s->native_map, s->a);
}

static void simulation_save_state(struct simulation *self)
{
        protein_print_atoms(self->protein, self->configurations);
        fflush(self->configurations);

        fprintf(self->energies, "%e\n", self->U);
        fflush(self->energies);
}

void simulation_do_iteration(struct simulation *self)
{
        assert(self != NULL);

        ++self->total;

        struct protein *candidate, *current = self->protein;

        candidate = protein_dup(current);
        protein_do_movements(candidate, self->rng);

        const double U1 = compute_potential_energy(current, self);
        const double U2 = compute_potential_energy(candidate, self);
        const double DU = U2 - U1;

        struct protein *chosen;

        if (DU < 0.0)
                chosen = candidate;
        else {
                const double T = 25.5;   /* XXX Change the temperature. */
                const double r = gsl_rng_uniform(self->rng);
                const double p = exp(-1.0/(GSL_CONST_MKSA_BOLTZMANN*T)*DU);

                chosen = r < p ? candidate : current;
        }

        if (chosen == candidate) {
                ++self->accepted;
                delete_protein(self->protein);
                self->protein = candidate;
                self->U = U2;
                simulation_save_state(self);
        } else {
                delete_protein(candidate);
        }
}
