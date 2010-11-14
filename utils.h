#ifndef UTILS_H
#define UTILS_H         1
/**
 * @file utils.h
 */

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


#define declare_stack_allocated_vector(name, n)                         \
        double stack_allocated_vector_array_##name[n];                  \
        gsl_vector_view stack_allocated_vector_view_##name =		\
                gsl_vector_view_array((double *) stack_allocated_vector_array_##name, n); \
        gsl_vector *name = &stack_allocated_vector_view_##name.vector

#define declare_stack_allocated_matrix(name, m, n)                      \
        double stack_allocated_matrix_array_##name[n];                  \
        gsl_matrix_view stack_allocated_matrix_view_##name =		\
                gsl_matrix_view_array((double *) stack_allocated_matrix_array_##name, m, n); \
        gsl_matrix *name = &stack_allocated_matrix_view_##name.matrix

#endif // !UTILS_H
