#ifndef PROTEIN_H
#define PROTEIN_H
/**
 * @file protein.h
 * @brief Interface for protein manipulation features.
 */

#include <stddef.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


/** Protein structure intended for coarse-grained models. */
struct protein {
        /** Number of alpha carbons. */
        size_t num_atoms;
        /** Array of vectors containing the position of each atom. */
        gsl_vector *atom[];
};


extern struct protein *new_protein(size_t num_atoms, const double *positions);
extern struct protein *new_protein_1pgb(void);
extern struct protein *new_protein_2gb1(void);
extern void delete_protein(struct protein *self);

extern int protein_print_atoms(const struct protein *self, FILE *stream);

extern double protein_distance(const struct protein *self, size_t i, size_t j);

extern void protein_plot(const struct protein *self);

extern void protein_move_end(struct protein *self, gsl_rng *rng);
#endif // !PROTEIN_H
