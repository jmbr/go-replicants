#ifndef REPLICAS_H
#define REPLICAS_H

struct replicas {
        gsl_rng *rng;

        /** Probability of exchanging replicas. */
        double alpha0;

        /** Number of replicas. */
        size_t num_replicas;

        /** Number of exchanges per pair of temperatures. */
        uint_fast32_t *exchanges;
        /** Number of attempted exchanges per temperature. */
        uint_fast32_t total;

        struct simulation *replica[];
};


extern struct replicas *new_replicas(const struct protein *p,
                                     gsl_rng *rng, double alpha0,
                                     size_t num_replicas,
                                     double temperatures[],
                                     const struct simulation_options *options);

extern void delete_replicas(struct replicas *self);

extern void replicas_first_iteration(struct replicas *self);
extern void replicas_next_iteration(struct replicas *self);

extern bool replicas_have_converged(const struct replicas *self);
#define replicas_have_not_converged(s)  (!replicas_have_converged(s))

extern uint_fast32_t replicas_total_exchanges(const struct replicas *self);

#endif // !REPLICAS_H
