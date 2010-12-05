#ifndef SIMULATION_H
#define SIMULATION_H

struct contact_map;

struct simulation {
        struct protein *protein;
        size_t next_atom;
        const struct contact_map *native_map;

        double a;

        double energy;

        double temperature;

        gsl_rng *rng;

        size_t accepted, total; /* Number of accepted moves and total number of moves. */

        FILE *U, *X;            /* Log files. */
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
