/**
 * @file simulation.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_const.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

#include "utils.h"
#include "protein.h"
#include "potential.h"
#include "contact-map.h"
#include "simulation.h"


const size_t simulation_save_step = 10000;

const char gnuplot_command_line[] = GNUPLOT_EXECUTABLE " -persist";



struct simulation *new_simulation(struct protein *p,
                                  const struct simulation_options *opts)
{
        if (p == NULL || opts == NULL || opts->d_max <= 0.0 || opts->a <= 0.0)
                return NULL;

        struct simulation *s;

        if ((s = calloc(1, sizeof(struct simulation))) == NULL)
                return NULL;

        s->protein = p;

        /* Initialize native contact map. */
        s->d_max = opts->d_max;
        if ((s->native_map = new_contact_map(s->protein, s->d_max)) == NULL)
                delete_simulation(s);

        s->a = opts->a;

        s->T = opts->T;

        /* Compute the potential energy for the protein. */
        s->orig_energy = s->energy = potential(s->protein, s->native_map, s->a);

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

static void simulation_save_state(const struct simulation *self)
{
        protein_print_atoms(self->protein, self->configurations);
        fflush(self->configurations);

        fprintf(self->energies, "%e\n", self->energy);
        fflush(self->energies);

        /* fputs("set terminal wxt noraise 1\n" */
        /*       "plot 'energies.dat' title 'Potential energy' with lines\n", */
        /*       self->gnuplot); */
}

void simulation_first_iteration(struct simulation *self)
{
        assert(self != NULL);

        ++self->total;

        protein_scramble(self->protein, self->rng);

        self->energy = compute_potential_energy(self->protein, self);

        simulation_save_state(self);
}

void simulation_next_iteration(struct simulation *self)
{
        assert(self != NULL);

        ++self->total;

        if (self->total % simulation_save_step == 0)
                simulation_save_state(self);

        struct protein *candidate, *current = self->protein;

        candidate = protein_dup(current);
        assert(candidate != NULL);
        protein_do_natural_movement(candidate, self->rng);

        const double U1 = compute_potential_energy(current, self);
        const double U2 = compute_potential_energy(candidate, self);
        const double DU = U2 - U1;

        dprintf("U1 == %g, U2 == %g, DU == %g\n", U1, U2, U2-U1);

        struct protein *chosen;

        if (DU <= 0.0)
                chosen = candidate;
        else {
                const double r = gsl_rng_uniform(self->rng);
                /* const double p = -1.0/(GSL_CONST_MKSA_BOLTZMANN*self->T)*DU; */
                const double p = -1.0/self->T*DU;

                chosen =  r < gsl_min(1.0, exp(p)) ? candidate : current;

                dprintf("r == %g, -1/T*DU == %g\n", r, p);
        }

        if (chosen == candidate) {
                ++self->accepted;
                delete_protein(self->protein);
                self->protein = candidate;
                self->energy = U2;
        } else {
                delete_protein(candidate);
        }
}

bool simulation_has_converged(const struct simulation *self)
{
        assert(self != NULL);

        return gsl_fcmp(self->orig_energy, self->energy, 1e-3) == 0;
}
