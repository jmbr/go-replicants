#ifndef PROTEIN_H
#define PROTEIN_H
/**
 * @file protein.h
 * @brief Interface for protein manipulation features.
 */

#include <stddef.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


/**
 * Protein structure intended for coarse-grained models.
 */
struct protein {
        size_t num_atoms;
        gsl_matrix *positions;
        gsl_matrix *native_contacts;
};


extern struct protein *new_protein(size_t num_atoms, const double *positions);
extern struct protein *new_protein_1pgb(void);
extern struct protein *new_protein_2gb1(void);

extern void delete_protein(struct protein *self);

extern size_t protein_num_atoms(const struct protein *self);

extern gsl_vector_view protein_get_atom(const struct protein *self, size_t index);

extern double protein_distance(const struct protein *self, size_t i, size_t j);

extern void protein_compute_contact_map(struct protein *self, double d_max, gsl_matrix *contact_map);

extern void protein_plot(const struct protein *self);

/* extern int protein_plot_contact_map(const struct protein *self); */
#endif // !PROTEIN_H
