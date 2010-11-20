#ifndef PROTEIN_H
#define PROTEIN_H
/**
 * @file protein.h
 */

#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


/** Coarse-grained protein structure. */
struct protein {
        /** Number of alpha carbons. */
        size_t num_atoms;
        /** Array of vectors containing the position of each atom. */
        gsl_vector *atom[];
};

enum protein_movements {
        PROTEIN_SPIKE_MOVE = 0,
        PROTEIN_SHIFT_MOVE,
        PROTEIN_PIVOT_MOVE,
        PROTEIN_END_MOVE_FIRST,
        PROTEIN_END_MOVE_LAST,
};


extern struct protein *new_protein(size_t num_atoms, const double *atom);
extern struct protein *new_protein_1pgb(void);
extern struct protein *new_protein_2gb1(void);
extern void delete_protein(struct protein *self);

extern struct protein *protein_dup(const struct protein *self);

extern double protein_distance(const struct protein *self, size_t i, size_t j);

extern void protein_plot(const struct protein *self, FILE *gnuplot,
                         const char *title_format, ...)
        __attribute__ ((format (printf, 3, 4)));

extern int protein_print_atoms(const struct protein *self, FILE *stream);

extern bool protein_is_overlapping(const struct protein *self);
#define protein_is_not_overlapping(p)   !protein_is_overlapping(p)

extern void protein_do_movement(struct protein *self, gsl_rng *rng,
                                enum protein_movements m, size_t k, bool undo);
extern void protein_do_natural_movement(struct protein *self, gsl_rng *rng);

extern void protein_scramble(struct protein *self, gsl_rng *rng);
#endif // !PROTEIN_H
