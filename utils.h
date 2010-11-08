#ifndef UTILS_H
#define UTILS_H         1
/**
 * @file utils.h
 */

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


extern void cross_product(const gsl_vector *u, const gsl_vector *v, gsl_vector *result);

extern double triple_scalar_product(const gsl_vector *u, const gsl_vector *v, const gsl_vector *w);

extern int print_matrix(FILE *stream, const gsl_matrix *matrix);

#endif // !UTILS_H
