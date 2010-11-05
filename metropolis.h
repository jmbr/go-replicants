#ifndef METROPOLIS_H
#define METROPOLIS_H
/*
 * metropolis.h -- Interface to the Metropolis-Hastings algorithm.
 */

struct metropolis;


typedef void *metropolis_state;

typedef metropolis_state (*metropolis_generator)(metropolis_state);

typedef double (*metropolis_pi)(metropolis_state);


extern struct metropolis *new_metropolis(unsigned int seed, metropolis_generator generator, metropolis_pi pi);

extern void delete_metropolis(struct metropolis *self);

extern double metropolis_random(const struct metropolis *self);

extern metropolis_state metropolis_iteration(const struct metropolis *self, metropolis_state current_state);
#endif // !METROPOLIS_H
