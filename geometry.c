#include "molecular-simulator.h"


static void make_random_unit_quaternion(gsl_vector *Q, gsl_rng *rng);
static void make_rotation_matrix_from_unit_quaternion(const gsl_vector *Q,
                                                      double R[3][3]);

void cross_product(const gsl_vector *u,
                   const gsl_vector *v,
                   gsl_vector *product)
{
        const double p1 = gsl_vector_get(u, 1)*gsl_vector_get(v, 2)
                - gsl_vector_get(u, 2)*gsl_vector_get(v, 1);

        const double p2 = gsl_vector_get(u, 2)*gsl_vector_get(v, 0)
                - gsl_vector_get(u, 0)*gsl_vector_get(v, 2);

        const double p3 = gsl_vector_get(u, 0)*gsl_vector_get(v, 1)
                - gsl_vector_get(u, 1)*gsl_vector_get(v, 0);

        gsl_vector_set(product, 0, p1);
        gsl_vector_set(product, 1, p2);
        gsl_vector_set(product, 2, p3);
}

double triple_scalar_product(const gsl_vector *u,
                             const gsl_vector *v,
                             const gsl_vector *w)
{
        double space[3], result;
        gsl_block b = block_init(3, space);
        gsl_vector vxw = vector_init(b, 0, 3, 1);

        cross_product(v, w, &vxw);
        gsl_blas_ddot(u, &vxw, &result);

        return result;
}

void make_random_rotation_matrix(double R[3][3], gsl_rng *rng)
{
        declare_stack_allocated_vector(Q, 4);
        make_random_unit_quaternion(Q, rng);

        make_rotation_matrix_from_unit_quaternion(Q, R);
}

void make_random_unit_quaternion(gsl_vector *Q, gsl_rng *rng)
{
        /* This comes from Graphics Gems III p. 129. */
        for (size_t i = 0; i < 4; i++)
                gsl_vector_set(Q, i, gsl_ran_ugaussian(rng));
        vector_normalize(Q);
}

void make_rotation_matrix_from_unit_quaternion(const gsl_vector *Q,
                                               double R[3][3])
{
        R[0][0] = gsl_pow_2(gsl_vector_get(Q, 0))
                + gsl_pow_2(gsl_vector_get(Q, 1))
                - gsl_pow_2(gsl_vector_get(Q, 2))
                - gsl_pow_2(gsl_vector_get(Q, 3));
        R[0][1] = 2.0*(gsl_vector_get(Q, 1)*gsl_vector_get(Q, 2)
                       + gsl_vector_get(Q, 0)*gsl_vector_get(Q, 3));
        R[0][2] = 2.0*(gsl_vector_get(Q, 1)*gsl_vector_get(Q, 3)
                       - gsl_vector_get(Q, 0)*gsl_vector_get(Q, 2));

        R[1][0] = 2.0*(gsl_vector_get(Q, 1)*gsl_vector_get(Q, 2)
                       - gsl_vector_get(Q, 0)*gsl_vector_get(Q, 3));
        R[1][1] = gsl_pow_2(gsl_vector_get(Q, 0))
                - gsl_pow_2(gsl_vector_get(Q, 1))
                + gsl_pow_2(gsl_vector_get(Q, 2))
                - gsl_pow_2(gsl_vector_get(Q, 3));
        R[1][2] = 2.0*(gsl_vector_get(Q, 2)*gsl_vector_get(Q, 3)
                       + gsl_vector_get(Q, 0)*gsl_vector_get(Q, 1));

        R[2][0] = 2.0*(gsl_vector_get(Q, 1)*gsl_vector_get(Q, 3)
                       + gsl_vector_get(Q, 0)*gsl_vector_get(Q, 2));
        R[2][1] = 2.0*(gsl_vector_get(Q, 2)*gsl_vector_get(Q, 3)
                       - gsl_vector_get(Q, 0)*gsl_vector_get(Q, 1));
        R[2][2] = gsl_pow_2(gsl_vector_get(Q, 0))
                - gsl_pow_2(gsl_vector_get(Q, 1))
                - gsl_pow_2(gsl_vector_get(Q, 2))
                + gsl_pow_2(gsl_vector_get(Q, 3));
}

/** Rotates point b around point a using the rotation matrix R. */
void rotate(bool transpose, const gsl_matrix *R, const gsl_vector *a, gsl_vector *b)
{
        declare_stack_allocated_vector(v, 3);
        gsl_vector_memcpy(v, b);
        gsl_vector_sub(v, a);

        /* Rotate end vector. */
        declare_stack_allocated_vector(w, 3);
        gsl_blas_dgemv(transpose == false ? CblasNoTrans : CblasTrans, 1.0, R, v, 0.0, w);

        /* Update position. */
        gsl_vector_memcpy(b, a);
        gsl_vector_add(b, w);

        assert(gsl_fcmp(gsl_blas_dnrm2(v), gsl_blas_dnrm2(w), 1e-15) == 0);
}

int print_matrix(FILE *stream, const gsl_matrix *matrix)
{
        int status, n = 0;

        for (size_t i = 0; i < matrix->size1; i++) {
                for (size_t j = 0; j < matrix->size2; j++) {
                        status = fprintf(stream, "%g ",
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
                status = fprintf(stream, "%g ", gsl_vector_get(vector, i));
                if (status < 0)
                        return -1;
                n += status;
        }
        if ((status = fprintf(stream, "\n")) < 0)
                return -1;
        n += status;

        return n;
}
