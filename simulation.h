#ifndef SIMULATION_H
#define SIMULATION_H

/**
 * Simulation structure.  Provides context for the simulation
 * procedure.
 */
struct simulation {
        /** Molecule being simulated. */
        struct protein *protein;

        /** Potential energy. */
        double orig_energy, energy;

        /** Temperature. */
        double T;

        /** Tolerance for distances between amino acids.  This is used
         * for the potential energy calculation and its value should
         * be less than 1 Angstrom. */
        double a;

        /** Maximum distance between native contacts. */
        double d_max;

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

        /** Index of the next atom from which to perform a
         * displacement. */
        size_t next_atom;
};

struct simulation_options {
        double d_max, a;
};


extern struct simulation *new_simulation(const struct protein *p,
                                         gsl_rng *rng,
                                         double temperature,
                                         const struct simulation_options *opts);
extern void delete_simulation(struct simulation *self);

extern void simulation_first_iteration(struct simulation *self);
extern void simulation_next_iteration(struct simulation *self);

extern bool simulation_has_converged(const struct simulation *self);
#define simulation_has_not_converged(s)         !simulation_has_converged(s)

extern double simulation_get_acceptance_ratio(const struct simulation *self);

extern void simulation_print_info(const struct simulation *self, FILE *stream);

#endif // !SIMULATION_H
