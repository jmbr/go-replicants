#ifndef PROTEIN_H
#define PROTEIN_H

#include "movements.h"


/** Coarse-grained protein structure. */
struct protein {
        /** Number of alpha carbons. */
        size_t num_atoms;
        /** Array of vectors containing the position of each atom. */
        gsl_vector *atom[];
};


/* Allocation and deallocation of protein data structures. */
extern struct protein *new_protein(size_t num_atoms, const double *atom);
extern struct protein *new_protein_1pgb(void);
extern struct protein *new_protein_2gb1(void);
extern void delete_protein(struct protein *self);
extern struct protein *protein_dup(const struct protein *self);

/* Input/Output functions. */
extern struct protein *protein_read_xyz_file(const char *name);
extern struct protein *protein_read_xyz(FILE *stream);
extern struct protein *protein_read_latest_xyz(FILE *stream);
extern int protein_write_xyz_file(const struct protein *self, const char *name);
extern int protein_write_xyz(const struct protein *self, FILE *stream);

extern int protein_print_atoms(const struct protein *self, FILE *stream);

extern void protein_plot(const struct protein *self, FILE *gnuplot,
                         bool draw_labels, const char *title_format, ...)
        __attribute__ ((format (printf, 4, 5)));

/* Geometry functions. */
extern double protein_signum(const struct protein *self, size_t i, size_t j);
extern double protein_distance(const struct protein *self, size_t i, size_t j);

#endif // !PROTEIN_H
