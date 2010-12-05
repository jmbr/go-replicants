#ifndef MOVEMENTS_H
#define MOVEMENTS_H

struct protein;

enum protein_movements {
        PROTEIN_SPIKE_MOVE = 0,
        PROTEIN_SHIFT_MOVE,
        PROTEIN_PIVOT_MOVE,
        PROTEIN_END_MOVE_FIRST,
        PROTEIN_END_MOVE_LAST,
};


extern bool protein_is_overlapping(const struct protein *self);
#define protein_is_not_overlapping(p)   !protein_is_overlapping(p)

extern bool protein_do_movement(struct protein *self, gsl_rng *rng,
                                enum protein_movements m, size_t k);
extern bool protein_do_natural_movement(struct protein *self, gsl_rng *rng, size_t k);

extern void protein_scramble(struct protein *self, gsl_rng *rng);

#endif // !MOVEMENTS_H
