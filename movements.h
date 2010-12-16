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


extern bool protein_is_overlapping(const struct protein *self, size_t start, size_t end);
#define protein_is_not_overlapping(p, s, e)   !protein_is_overlapping(p, s, e)

extern bool protein_do_movement(struct protein *self, gsl_rng *rng,
                                enum protein_movements m, size_t k);

extern bool protein_do_natural_movement(struct protein *self, gsl_rng *rng, size_t k);

extern bool protein_do_shift_move(struct protein *self, gsl_rng *rng, size_t k);
extern bool protein_do_spike_move(struct protein *self, gsl_rng *rng, size_t k);
extern bool protein_do_pivot_move(struct protein *self, gsl_rng *rng, size_t k);
extern bool protein_do_end_move_first(struct protein *self, gsl_rng *rng);
extern bool protein_do_end_move_last(struct protein *self, gsl_rng *rng);

extern void protein_scramble(struct protein *self, gsl_rng *rng);

#endif // !MOVEMENTS_H
