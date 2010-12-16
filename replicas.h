#ifndef REPLICAS_H
#define REPLICAS_H

struct contact_map;

/** Replica data structure.  This provides context for the replica
 * exchange simulation. */
struct replicas {
        gsl_rng *rng;		        /**< Random number generator. */
        struct protein *protein;	/**< Protein to be simulated. */
        struct contact_map *native_map; /**< Native contacts. */
        double a;                       /**< Tolerance for distances between amino acids. */
        size_t num_replicas;            /**< Number of replicas. */
        size_t *exchanges;              /**< Number of exchanges per pair of replicas. */
        size_t *total;                  /**< Number of attempted exchanges per pair of replicas. */
        struct simulation *replica[];   /**< Array of replicas. */
};

/** Auxiliary options for new_replicas. */
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

extern size_t replicas_total_exchanges(const struct replicas *self);
extern void replicas_get_exchange_ratios(const struct replicas *self,
                                         double ratios[]);
extern void replicas_print_info(const struct replicas *self, FILE *stream);

#endif // !REPLICAS_H
