#ifndef METROPOLIS_H
#define METROPOLIS_H
/**
 * @file metropolis.h
 * @brief Interface to the Metropolis-Hastings algorithm.
 */

struct metropolis;


typedef void *metropolis_state;

typedef metropolis_state (*metropolis_generator)(metropolis_state);

typedef double (*metropolis_pi)(metropolis_state);


extern struct metropolis *new_metropolis(unsigned int seed, metropolis_generator generator, metropolis_pi pi);

extern void delete_metropolis(struct metropolis *self);

extern double metropolis_random(const struct metropolis *self);

extern metropolis_state metropolis_iteration(struct metropolis *self, metropolis_state current_state);

extern double metropolis_get_acceptance_ratio(const struct metropolis *self);
#endif // !METROPOLIS_H
