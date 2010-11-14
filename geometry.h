#ifndef GEOMETRY_H
#define GEOMETRY_H	1
/**
 * @file geometry.h
 */

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>


extern void cross_product(const gsl_vector *u, const gsl_vector *v,
                          gsl_vector *result);

extern double triple_scalar_product(const gsl_vector *u,
                                    const gsl_vector *v,
                                    const gsl_vector *w);

extern int print_matrix(FILE *stream, const gsl_matrix *matrix);

extern int print_vector(FILE *stream, const gsl_vector *vector);

static inline void vector_normalize(gsl_vector *v)
{
        gsl_vector_scale(v, 1.0/gsl_blas_dnrm2(v));
}
#endif // !GEOMETRY_H
