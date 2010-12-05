#ifndef REPLICAS_H
#define REPLICAS_H

struct contact_map;

struct replicas {
        gsl_rng *rng;

        struct protein *protein;

        /** Tolerance for distances between amino acids.  This is used
         * for the potential energy calculation and its value should
         * be less than 1 Angstrom. */
        double a;

        struct contact_map *native_map;
        double orig_energy;

        /** Probability of exchanging replicas. */
        double alpha0;

        /** Number of replicas. */
        size_t num_replicas;

        /** Number of exchanges per pair of temperatures and number of
         * attempted exchanges per pair of temperature. */
        size_t *exchanges, *total;

        struct simulation *replica[];
};

struct simulation_options {
        gsl_rng *rng;
        double d_max, a;
        size_t num_replicas;
        double *temperatures;
};


extern struct replicas *new_replicas(struct protein *protein,
                                     const struct simulation_options *options);

extern void delete_replicas(struct replicas *self);

extern void replicas_thermalize(struct replicas *self, size_t num_iters);
extern void replicas_resume(struct replicas *self,
                            const struct protein *config[]);
extern void replicas_first_iteration(struct replicas *self);
extern void replicas_next_iteration(struct replicas *self);

extern bool replicas_have_converged(const struct replicas *self);
#define replicas_have_not_converged(s)  (!replicas_have_converged(s))

extern size_t replicas_total_exchanges(const struct replicas *self);
extern void replicas_get_exchange_ratios(const struct replicas *self,
                                         double ratios[]);
extern void replicas_print_info(const struct replicas *self, FILE *stream);

#endif // !REPLICAS_H
