#ifndef GEOMETRY_H
#define GEOMETRY_H	1
/**
 * @file geometry.h
 */

#include <stdbool.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>


extern void cross_product(const gsl_vector *u, const gsl_vector *v,
                          gsl_vector *result);

extern double triple_scalar_product(const gsl_vector *u,
                                    const gsl_vector *v,
                                    const gsl_vector *w);

extern void make_random_rotation_matrix(double R[3][3], gsl_rng *rng);

extern void rotate(bool transpose, const gsl_matrix *R, const gsl_vector *a, gsl_vector *b);

extern int print_matrix(FILE *stream, const gsl_matrix *matrix);

extern int print_vector(FILE *stream, const gsl_vector *vector);
#if DEBUG_LEVEL > 0
# define dprint_vector(v) print_vector(stderr, v)
#else
# define dprint_vector(v)
#endif // DEBUG_LEVEL > 0

static inline void vector_normalize(gsl_vector *v)
{
        gsl_vector_scale(v, 1.0/gsl_blas_dnrm2(v));
}
#endif // !GEOMETRY_H
