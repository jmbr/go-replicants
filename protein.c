/**
 * @file protein.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdarg.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "utils.h"
#include "geometry.h"
#include "protein.h"



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
                m->atom[i] = gsl_vector_alloc(3);
                /* XXX GSL takes care of the error checking. */
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



static double signum(const struct protein *self, size_t i, size_t j)
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



void protein_plot(const struct protein *self, FILE *gnuplot, const char *title_format, ...)
{
        assert(self != NULL);
        assert(gnuplot != NULL);

        va_list ap;

        fprintf(gnuplot, "set title '");
        va_start(ap, title_format);
        vfprintf(gnuplot, title_format, ap);
        va_end(ap);
        fprintf(gnuplot, "'\n");

        fprintf(gnuplot, "set view equal xyz\n"
                         "set linetype 1 linecolor palette z linewidth 5\n"
                         "unset tics\n"
                         "unset border\n"
                         "unset colorbox\n"
                         "splot '-' title '' with lines\n");

        protein_print_atoms(self, gnuplot);

        fprintf(gnuplot, "e\n");

        fflush(gnuplot);
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



bool protein_is_overlapping(const struct protein *self)
{
        const double bead_diameter = 1.1*3.8;

        for (size_t i = 0; i < self->num_atoms; i++)
                for (size_t j = i+2; j < self->num_atoms; j++)
                        if (fabs(protein_distance(self, i, j)) < bead_diameter) {
                                dprintf("overlap of %g between atoms %u and %u\n",
                                        fabs(protein_distance(self, i, j)), i, j);
                                return true;
                        }


        return false;
}



void protein_do_movements(struct protein *self, gsl_rng *rng)
{
        /* TODO: Prove that no redundant movements are being done. */
        protein_do_end_move_first(self, rng);
        protein_do_spike_move(self, rng, 0);
        for (size_t i = 1; i < self->num_atoms - 2; i++) {
                protein_do_spike_move(self, rng, i);
                protein_do_shift_move(self, rng, i);
                protein_do_pivot_move(self, rng, i);
        }
        protein_do_end_move_last(self, rng);
}



void protein_do_end_move_first(struct protein *self, gsl_rng *rng)
{
        assert(self != NULL);
        assert(self->num_atoms >= 2);
        assert(rng != NULL);

        dprintf("moving first atom.\n");
        dprintf("before: atom(0) == "); print_vector(stderr, self->atom[0]);

        double R[3][3];
        make_random_rotation_matrix(R, rng);
        gsl_matrix_view RV = gsl_matrix_view_array((double *) R, 3, 3);
        rotate(false, &RV.matrix, self->atom[1], self->atom[0]);
        dprintf("after: atom(0) == "); print_vector(stderr, self->atom[0]);

        if (protein_is_overlapping(self)) {
                /* Let's undo the movement. */
                rotate(true, &RV.matrix, self->atom[1], self->atom[0]);
                dprintf("after undo: atom(0) == "); print_vector(stderr, self->atom[0]);
        }
}

void protein_do_end_move_last(struct protein *self, gsl_rng *rng)
{
        assert(self != NULL);
        assert(self->num_atoms >= 2);
        assert(rng != NULL);

        const size_t num_atoms = self->num_atoms;

        dprintf("moving last atom.\n");
        dprintf("before: atom(%d) == ", num_atoms-1); print_vector(stderr, self->atom[num_atoms - 1]);

        double R[3][3];
        make_random_rotation_matrix(R, rng);
        gsl_matrix_view RV = gsl_matrix_view_array((double *) R, 3, 3);
        rotate(false, &RV.matrix, self->atom[self->num_atoms - 2],
                           self->atom[self->num_atoms - 1]);
        dprintf("after: atom(%d) == ", num_atoms-1); print_vector(stderr, self->atom[num_atoms - 1]);

        if (protein_is_overlapping(self)) {
                /* Let's undo the movement. */
                rotate(true, &RV.matrix,
                       self->atom[self->num_atoms - 2],
                       self->atom[self->num_atoms - 1]);
                dprintf("after undo: atom(%d) == ", num_atoms-1); print_vector(stderr, self->atom[num_atoms-1]);
        }
}



void protein_do_shift_move(struct protein *self, gsl_rng *rng, size_t k)
{
        assert(self != NULL);
        assert(self->num_atoms >= 3);
        assert(rng != NULL);

        assert(k <= self->num_atoms - 3);

        /* const size_t k = gsl_rng_uniform_int(rng, self->num_atoms - 3); */

        dprintf("shifting move at atom %u\n", k);
        dprintf("before: atom(%u) == ", self->num_atoms-1); print_vector(stderr, self->atom[self->num_atoms-1]);
        dprintf("before: atom(%u) == ", k+1); print_vector(stderr, self->atom[k+1]);

        double R[3][3];
        make_random_rotation_matrix(R, rng);
        gsl_matrix_const_view RV = gsl_matrix_const_view_array((double *) R, 3, 3);

        /*
         * 1. Take a consecutive pair of atoms: a(k) and a(k+1).
         * 2. Apply a random rotation to the vector joining them.
         * 3. Translate the atoms a(k+2) ... a(num_atoms-1) so that
         * the connectivity of the chain is maintained.
         */
        gsl_vector *bak = gsl_vector_alloc(3);
        gsl_vector_memcpy(bak, self->atom[k+1]);

        /* t = atom(k+1) - R atom(k+1) */
        declare_stack_allocated_vector(t, 3);
        gsl_vector_memcpy(t, self->atom[k+1]);
        rotate(false, &RV.matrix, self->atom[k], self->atom[k+1]);
        gsl_vector_sub(t, self->atom[k+1]);

        for (size_t i = k+2; i < self->num_atoms; i++)
                gsl_vector_sub(self->atom[i], t);

        dprintf("after: atom(%u) == ", self->num_atoms-1); print_vector(stderr, self->atom[self->num_atoms-1]);
        dprintf("after: atom(%u) == ", k+1); print_vector(stderr, self->atom[k+1]);

        if (protein_is_overlapping(self)) { /* Undo movement. */
                dprintf("undoing shift movement.\n");

                for (size_t i = k+2; i < self->num_atoms; i++)
                        gsl_vector_add(self->atom[i], t);

                dprintf("after undo: atom(%u) == ", self->num_atoms-1); print_vector(stderr, self->atom[self->num_atoms-1]);

                gsl_vector_memcpy(self->atom[k+1], bak);

                dprintf("after undo: atom(%u) == ", k+1); print_vector(stderr, self->atom[k+1]);
        }

        gsl_vector_free(bak);
}



static void do_spike_move(struct protein *self, size_t i, double theta)
{
        dprintf("performing spike move of %g radians on atom %u.\n", theta, i);

        gsl_vector *bak = gsl_vector_alloc(3);
        gsl_vector_memcpy(bak, self->atom[i+1]);

        dprintf("before: atom(%u) == ", i+1); print_vector(stderr, self->atom[i+1]);

        /* v = p3-p1 */
        gsl_vector *v = gsl_vector_alloc(3);
        gsl_vector_memcpy(v, self->atom[i+2]);
        gsl_vector_sub(v, self->atom[i]);
        vector_normalize(v);

        /* a = p2-p1 */
        gsl_vector *a = gsl_vector_alloc(3);
        gsl_vector_memcpy(a, self->atom[i+1]);
        gsl_vector_sub(a, self->atom[i]);

        /* w = <a, v> v */
        gsl_vector *w = gsl_vector_alloc(3);
        gsl_vector_memcpy(w, v);
        double len;
        gsl_blas_ddot(a, v, &len);
        gsl_vector_scale(w, len);

        /* t = p1 + w */
        gsl_vector *t = gsl_vector_alloc(3);
        gsl_vector_memcpy(t, w);
        gsl_vector_add(t, self->atom[i]);

        /* q = atom[i+1] - t */
        gsl_vector *q = gsl_vector_alloc(3);
        gsl_vector_memcpy(q, self->atom[i+1]);
        gsl_vector_sub(q, t);

        gsl_matrix *G = gsl_matrix_alloc(3, 3);
        gsl_vector_view u0 = gsl_matrix_column(G, 0);
        gsl_vector_view u1 = gsl_matrix_column(G, 1);
        gsl_vector_view u2 = gsl_matrix_column(G, 2);

        gsl_vector_memcpy(&u0.vector, q);
        vector_normalize(&u0.vector);

        gsl_vector_memcpy(&u2.vector, w);
        vector_normalize(&u2.vector);

        cross_product(&u0.vector, &u2.vector, &u1.vector);
        if (gsl_fcmp(gsl_blas_dnrm2(&u1.vector), 1.0, GSL_DBL_EPSILON) != 0) {
                dprintf("|u0 x u2| == "); print_vector(stderr, &u1.vector);
                dprintf("ignoring collinear configuration:\n");
                dprintf("atom(%u) == ", i);   print_vector(stderr, self->atom[i]);
                dprintf("atom(%u) == ", i+1); print_vector(stderr, self->atom[i+1]);
                dprintf("atom(%u) == ", i+2); print_vector(stderr, self->atom[i+2]);
                fflush(stderr);

                assert(0);
                goto cleanup;
        }

        double R[3][3] = {{cos(theta), -sin(theta), 0.0},
                          {sin(theta),  cos(theta), 0.0},
                          {       0.0,         0.0, 1.0}};
        gsl_matrix_const_view RV = gsl_matrix_const_view_array((double *) R, 3, 3);

        gsl_matrix *A = gsl_matrix_alloc(3, 3);
        gsl_matrix *B = gsl_matrix_alloc(3, 3);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &RV.matrix, G, 0.0, A);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, G, A, 0.0, B);

        gsl_vector *z = gsl_vector_alloc(3);
        gsl_blas_dgemv(CblasNoTrans, 1.0, B, q, 0.0, z);

        gsl_vector_add(z, t);
        gsl_vector_memcpy(self->atom[i+1], z);

        dprintf("after: atom(%d) == ", i+1); print_vector(stderr, self->atom[i+1]);

        /* We are done if the conformation is correct. */
        if (protein_is_overlapping(self)) {
                dprintf("undoing invalid conformation.\n");

                gsl_vector_memcpy(self->atom[i+1], bak);

                dprintf("after undo: atom(%d) == ", i+1); print_vector(stderr, self->atom[i+1]);
        }

cleanup:                        /* XXX This is a mess. */
        gsl_vector_free(bak);
        gsl_vector_free(v);
        gsl_vector_free(a);
        gsl_vector_free(w);
        gsl_vector_free(t);
        gsl_vector_free(q);
        gsl_vector_free(z);
        gsl_matrix_free(G);
        gsl_matrix_free(A);
        gsl_matrix_free(B);
}

void protein_do_spike_move(struct protein *self, gsl_rng *rng, size_t k)
{
        assert(self != NULL);
        assert(self->num_atoms >= 3);
        assert(rng != NULL);

        if (k > self->num_atoms - 3) {
        }

        const double theta = 2*M_PI*gsl_rng_uniform_pos(rng);
        /* const size_t k = gsl_rng_uniform_int(rng, self->num_atoms - 3); */

        do_spike_move(self, k, theta);
}



void protein_do_pivot_move(struct protein *self, gsl_rng *rng, size_t k)
{
        assert(self != NULL);
        assert(rng != NULL);

        /* const size_t k = gsl_rng_uniform_int(rng, self->num_atoms - 1); */
        const double theta = 2*M_PI*gsl_rng_uniform_pos(rng);

        dprintf("pivoting move at atom %u.\n", k);

        dprintf("before: atom(%d) == ", k+1);
        print_vector(stderr, self->atom[k+1]);

        double RR[3][3] = {{cos(theta), -sin(theta), 0.0},
                           {sin(theta),  cos(theta), 0.0},
                           {       0.0,         0.0, 1.0}};
        gsl_matrix_const_view RV = gsl_matrix_const_view_array((double *) RR, 3, 3);
        declare_stack_allocated_vector(y, 3);
        for (size_t i = k+1; i < self->num_atoms; i++) {
                gsl_vector_memcpy(y, self->atom[k]);
                gsl_vector_sub(self->atom[i], y);
                gsl_blas_dgemv(CblasNoTrans, 1.0, &RV.matrix, self->atom[i],
                               1.0, y);
                gsl_vector_memcpy(self->atom[i], y);
        }

        dprintf("after: atom(%d) == ", k+1);
        print_vector(stderr, self->atom[k+1]);

        /* We are done if the conformation is correct. */
        if (protein_is_not_overlapping(self))
                return;

        dprintf("undoing invalid conformation.\n");

        /* Undo movement if an invalid conformation is reached. */
        for (size_t i = k+1; i < self->num_atoms; i++) {
                gsl_vector_memcpy(y, self->atom[k]);
                gsl_vector_sub(self->atom[i], y);
                gsl_blas_dgemv(CblasTrans, 1.0, &RV.matrix, self->atom[i],
                               1.0, y);
                gsl_vector_memcpy(self->atom[i], y);
        }

        dprintf("after undo: atom(%d) == ", k+1);
        print_vector(stderr, self->atom[k+1]);
}
