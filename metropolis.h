#ifndef METROPOLIS_H
#define METROPOLIS_H
/*
 * metropolis.h -- Interface to the Metropolis-Hastings algorithm.
 */

struct metropolis;


typedef void *state;

typedef state (*generator_fn)(state);

typedef double (*pi_fn)(state);


extern struct metropolis *new_metropolis(unsigned int seed, generator_fn generator, pi_fn pi);

extern void delete_metropolis(struct metropolis *self);

extern double metropolis_random(const struct metropolis *self);

extern state metropolis_iteration(const struct metropolis *self, state current_state);
#endif // !METROPOLIS_H
