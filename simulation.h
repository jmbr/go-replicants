#ifndef SIMULATION_H
#define SIMULATION_H            1
/**
 * @file simulation.h
 */

#include <stdio.h>
#include <stdint.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>


struct protein;

struct contact_map;

/**
 * Simulation structure.  Provides context for the simulation
 * procedure.
 */
struct simulation {
        /** Molecule being simulated. */
        struct protein *protein;

        /** Potential energy. */
        double U;

        /** Maximum distance between native contacts. */
        double d_max;

        /** Contact map for the native configuration. */
        struct contact_map *native_map;

        /** Random number generator. */
        gsl_rng *rng;

        /** Tolerance for distances between amino acids.  This is used
         * for the potential energy calculation and its value should
         * be less than 1 Angstrom. */
        double a;

        /** Number of accepted moves and number of total moves. */
        uint_fast32_t accepted, total;

        /** Instance of Gnuplot. */
        FILE *gnuplot;

        /** Files containing the results. */
        FILE *configurations, *energies;
};


extern struct simulation *new_simulation(struct protein *p,
                                         double d_max, double a);
extern void delete_simulation(struct simulation *self);

extern void simulation_do_iteration(struct simulation *self);

extern double simulation_get_acceptance_ratio(const struct simulation *self);
#endif // !SIMULATION_H
