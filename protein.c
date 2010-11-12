/**
 * @file protein.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "utils.h"
#include "protein.h"


static const char gnuplot_command_line[] = GNUPLOT_EXECUTABLE " --persist";


static double signum(const struct protein *self, size_t i, size_t j);
static void make_uniformly_distributed_random_rotation_matrix(double R[3][3], gsl_rng *rng);
static void make_uniformly_distributed_random_unit_quaternion(double Q[4], gsl_rng *rng);
static void make_rotation_matrix_from_unit_quaternion(const double Q[4], double R[3][3]);



/////////////////////////////////////////////////////////////////////////////
// Allocation and deallocation of protein structures.
/////////////////////////////////////////////////////////////////////////////

struct protein *new_protein(size_t num_atoms, const double *atom)
{
        struct protein *m;

        if (num_atoms == 0 || atom == NULL)
                return NULL;

        m = malloc(sizeof(struct protein) + num_atoms*sizeof(gsl_vector *));
        if (m == NULL)
                return NULL;

        m->num_atoms = num_atoms;

        for (size_t i = 0; i < m->num_atoms; i++) {
                m->atom[i] = gsl_vector_alloc(3); /* XXX GSL takes care of the error checking. */
                for (size_t j = 0; j < 3; j++)
                        gsl_vector_set(m->atom[i], j, atom[3*i + j]);
        }

        return m;
}

#include "sample-proteins.c"

void delete_protein(struct protein *self)
{
        assert(self != NULL);
        assert(self->atom != NULL);

        for (size_t i = 0; i < self->num_atoms; i++)
                gsl_vector_free(self->atom[i]);
        free(self);
}



double protein_distance(const struct protein *self,
                        size_t i, size_t j)
{
        assert(self != NULL);
        assert(i < self->num_atoms);
        assert(j < self->num_atoms);

        double data[] = {0, 0, 0};
        gsl_vector_view v = gsl_vector_view_array((double *) &data, 3);
        gsl_vector_memcpy(&v.vector, self->atom[j]);
        gsl_vector_sub(&v.vector, self->atom[i]);

        return signum(self, i, j)*gsl_blas_dnrm2(&v.vector);
}

double signum(const struct protein *self, size_t i, size_t j)
{
        double s = 1.0;

        if (abs((int) (i - j)) == 3) {
                gsl_vector *v0 = self->atom[i];
                gsl_vector *v1 = self->atom[i+1];
                gsl_vector *v2 = self->atom[i+2];
                gsl_vector *v3 = self->atom[i+3];

                double uu[3], vv[3], ww[3];
                gsl_vector_view u = gsl_vector_view_array((double *) &uu, 3);
                gsl_vector_view v = gsl_vector_view_array((double *) &vv, 3);
                gsl_vector_view w = gsl_vector_view_array((double *) &ww, 3);

                gsl_vector_memcpy(&u.vector, v1);
                gsl_vector_memcpy(&v.vector, v2);
                gsl_vector_memcpy(&w.vector, v3);

                gsl_vector_sub(&u.vector, v0);
                gsl_vector_sub(&v.vector, v1);
                gsl_vector_sub(&w.vector, v2);

                double d = triple_scalar_product(&u.vector, &v.vector, &w.vector);

                s = d/fabs(d);
        }

        return s;
}



void protein_plot(const struct protein *self)
{
        assert(self != NULL);
        assert(self->atom != NULL);

        static FILE *g = NULL;

        if (g == NULL)
                g = popen(gnuplot_command_line, "w");
        /* if (g == NULL) */
        /*         return -1; */

        fprintf(g, "set title 'Alpha carbon structure'\n"
                   "set view equal\n"
                   "set linetype 1 linecolor palette z linewidth 5\n"
                   "unset tics\n"
                   "unset border\n"
                   "unset colorbox\n"
                   "splot '-' title '' with lines\n");

        protein_print_atoms(self, g);

        fprintf(g, "e\n");

        fflush(g);
        /* pclose(g); */
}

int protein_print_atoms(const struct protein *self, FILE *stream)
{
        int status, n = 0;

        for (size_t i = 0; i < self->num_atoms; i++) {
                status = fprintf(stream, "%e %e %e\n",
                                 gsl_vector_get(self->atom[i], 0),
                                 gsl_vector_get(self->atom[i], 1),
                                 gsl_vector_get(self->atom[i], 2));
                if (status < 0)
                        return -1;
                n += status;
        }

        return n;
}



/////////////////////////////////////////////////////////////////////////////
// Functions for performing changes in the protein's spatial structure.
/////////////////////////////////////////////////////////////////////////////

void protein_move_end(struct protein *self, gsl_rng *rng)
{
        assert(self != NULL);
        assert(rng != NULL);

        gsl_vector *a, *b;

        if (gsl_rng_get(rng) % 2UL == 0UL) {
                a = self->atom[self->num_atoms - 2];
                b = self->atom[self->num_atoms - 1];
        } else {
                a = self->atom[1];
                b = self->atom[0];
        }

        /* Fetch end vector. */
        double V[3];
        gsl_vector_view vv = gsl_vector_view_array(V, 3);
        gsl_vector_memcpy(&vv.vector, b);
        gsl_vector_sub(&vv.vector, a);

        /* Rotate end vector. */
        double R[3][3], W[3];
        make_uniformly_distributed_random_rotation_matrix(R, rng);
        gsl_matrix_const_view RV = gsl_matrix_const_view_array((double *) R, 3, 3);
        gsl_vector_view wv = gsl_vector_view_array(W, 3);
        gsl_blas_dgemv(CblasNoTrans, 1.0, &RV.matrix, &vv.vector, 0.0, &wv.vector);

        /* Update position. */
        gsl_vector_memcpy(b, a);
        gsl_vector_add(b, &wv.vector);

        assert(gsl_fcmp(gsl_blas_dnrm2(&vv.vector),
                        gsl_blas_dnrm2(&wv.vector), 1e-15) == 0);

        protein_plot(self);
}

void make_uniformly_distributed_random_rotation_matrix(double R[3][3], gsl_rng *rng)
{
        double Q[4];
        make_uniformly_distributed_random_unit_quaternion(Q, rng);
        make_rotation_matrix_from_unit_quaternion(Q, R);
}

void make_uniformly_distributed_random_unit_quaternion(double Q[4], gsl_rng *rng)
{
        /* This comes from Graphics Gems III p. 129. */
        for (size_t i = 0; i < 4; i++)
                Q[i] = gsl_ran_ugaussian(rng);
        gsl_vector_view QV = gsl_vector_view_array((double *) Q, 4);
        vector_normalize(&QV.vector);
}

void make_rotation_matrix_from_unit_quaternion(const double Q[4], double R[3][3])
{
        R[0][0] = gsl_pow_2(Q[0]) + gsl_pow_2(Q[1]) - gsl_pow_2(Q[2]) - gsl_pow_2(Q[3]);
        R[0][1] = 2.0*(Q[1]*Q[2] + Q[0]*Q[3]);
        R[0][2] = 2.0*(Q[1]*Q[3] - Q[0]*Q[2]);

        R[1][0] = 2.0*(Q[1]*Q[2] - Q[0]*Q[3]);
        R[1][1] = gsl_pow_2(Q[0]) - gsl_pow_2(Q[1]) + gsl_pow_2(Q[2]) - gsl_pow_2(Q[3]);
        R[1][2] = 2.0*(Q[2]*Q[3] + Q[0]*Q[1]);

        R[2][0] = 2.0*(Q[1]*Q[3] + Q[0]*Q[2]);
        R[2][1] = 2.0*(Q[2]*Q[3] - Q[0]*Q[1]);
        R[2][2] = gsl_pow_2(Q[0]) - gsl_pow_2(Q[1]) - gsl_pow_2(Q[2]) + gsl_pow_2(Q[3]);
}
