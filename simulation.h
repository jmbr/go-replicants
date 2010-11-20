#ifndef SIMULATION_H
#define SIMULATION_H            1
/**
 * @file simulation.h
 */

#include <stdio.h>
#include <stdbool.h>
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
        double orig_energy, energy;

        /** Maximum distance between native contacts. */
        double d_max;

        /** Tolerance for distances between amino acids.  This is used
         * for the potential energy calculation and its value should
         * be less than 1 Angstrom. */
        double a;

        /** Contact map for the native configuration. */
        struct contact_map *native_map;

        /** Random number generator. */
        gsl_rng *rng;

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

extern void simulation_first_iteration(struct simulation *self);
extern void simulation_next_iteration(struct simulation *self);

extern bool simulation_has_converged(const struct simulation *self);
#define simulation_has_not_converged(s)         !simulation_has_converged(s)

extern double simulation_get_acceptance_ratio(const struct simulation *self);
#endif // !SIMULATION_H
