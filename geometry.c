/**
 * @file geometry.h
 */

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "utils.h"
#include "geometry.h"


void cross_product(const gsl_vector *u,
                   const gsl_vector *v,
                   gsl_vector *product)
{
        double p1 = gsl_vector_get(u, 1)*gsl_vector_get(v, 2)
                - gsl_vector_get(u, 2)*gsl_vector_get(v, 1);

        double p2 = gsl_vector_get(u, 2)*gsl_vector_get(v, 0)
                - gsl_vector_get(u, 0)*gsl_vector_get(v, 2);

        double p3 = gsl_vector_get(u, 0)*gsl_vector_get(v, 1)
                - gsl_vector_get(u, 1)*gsl_vector_get(v, 0);

        gsl_vector_set(product, 0, p1);
        gsl_vector_set(product, 1, p2);
        gsl_vector_set(product, 2, p3);
}

double triple_scalar_product(const gsl_vector *u,
                             const gsl_vector *v,
                             const gsl_vector *w)
{
        double result;
        declare_stack_allocated_vector(vxw, 3);

        cross_product(v, w, vxw);
        gsl_blas_ddot(u, vxw, &result);

        return result;
}

int print_matrix(FILE *stream, const gsl_matrix *matrix)
{
        int status, n = 0;

        for (size_t i = 0; i < matrix->size1; i++) {
                for (size_t j = 0; j < matrix->size2; j++) {
                        status = fprintf(stream, "%e ",
                                         gsl_matrix_get(matrix, i, j));
                        if (status < 0)
                                return -1;
                        n += status;
                }

                if ((status = fprintf(stream, "\n")) < 0)
                        return -1;
                n += status;
        }

        return n;
}

int print_vector(FILE *stream, const gsl_vector *vector)
{
        int status, n = 0;

        for (size_t i = 0; i < vector->size; i++) {
                status = fprintf(stream, "%e ", gsl_vector_get(vector, i));
                if (status < 0)
                        return -1;
                n += status;
        }
        if ((status = fprintf(stream, "\n")) < 0)
                return -1;
        n += status;

        return n;
}
