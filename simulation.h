#ifndef SIMULATION_H
#define SIMULATION_H

struct contact_map;

/*** Individual replica.  This structure characterizes the simulation
 * process to be performed by an individual replica. */
struct simulation {
        struct protein *protein;                /**< Protein to be simulated. */
        size_t next_atom;			/**< Index of the next atom to be changed by a movement. */
        const struct contact_map *native_map;   /**< Native contacts. */
        double a;       		        /**< Tolerance parameter for the potential. */
        double energy;				/**< Current potential energy. */
        double temperature;                     /**< Temperature. */
        gsl_rng *rng;                           /**< Random number generator. */
        size_t accepted;                        /**< Number of accepted movements. */
        size_t total;                           /**< Number of attempted movements. */
        FILE *U;                                /**< Storage file containing energy values. */
        FILE *X;                                /**< Storage file containing spatial conformations. */
};


extern struct simulation *new_simulation(const struct contact_map *native_map,
                                         double a, double temperature,
                                         gsl_rng *rng);
extern void delete_simulation(struct simulation *self);

extern void simulation_first_iteration(struct simulation *self,
                                      const struct protein *protein, double energy);
extern void simulation_next_iteration(struct simulation *self);

extern double simulation_get_acceptance_ratio(const struct simulation *self);

extern void simulation_print_info(const struct simulation *self, FILE *stream);

#endif // !SIMULATION_H
