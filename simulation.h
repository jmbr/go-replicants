#ifndef SIMULATION_H
#define SIMULATION_H            1
/**
 * @file simulation.h
 */

#include <stdio.h>

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

        /** Instance of Gnuplot. */
        FILE *gnuplot;

        /* /\** Log file. *\/ */
        /* FILE *log; */
};

enum simulation_protein {
        SIMULATION_1PGB,
        SIMULATION_2GB1
};

extern struct simulation *new_simulation(enum simulation_protein p,
                                         double d_max, double a);

extern void delete_simulation(struct simulation *self);

/* extern double simulation_compute_potential(const simulation *self); */
#endif // !SIMULATION_H
