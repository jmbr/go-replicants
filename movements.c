#include "molecular-simulator.h"


bool protein_do_movement(struct protein *self, gsl_rng *rng,
                         enum protein_movements m, size_t k)
{
        bool status = false;

        switch (m) {
        case PROTEIN_SPIKE_MOVE:
                status = protein_do_spike_move(self, rng, k);
                break;
        case PROTEIN_SHIFT_MOVE:
                status = protein_do_shift_move(self, rng, k);
                break;
        case PROTEIN_PIVOT_MOVE:
                status = protein_do_pivot_move(self, rng, k);
                break;
        case PROTEIN_END_MOVE_FIRST:
                status = protein_do_end_move_first(self, rng);
                break;
        case PROTEIN_END_MOVE_LAST:
                status = protein_do_end_move_last(self, rng);
                break;
        }

        return status;
}

bool protein_do_natural_movement(struct protein *self, gsl_rng *rng, size_t k)
{
        /* dprintf("thread #%d is changing atom %u\n", omp_get_thread_num(), k); */

        if (k == 0) {
                return protein_do_end_move_first(self, rng);
        } else if (1 <= k && k <= self->num_atoms-2) {
                return protein_do_movement(self, rng, gsl_rng_uniform_int(rng, 2), k);
        } else {
                return protein_do_end_move_last(self, rng);
        }
}

bool protein_is_overlapping(const struct protein *self, size_t start, size_t end)
{
        const size_t N = self->num_atoms;
        const double bead_diameter = 1.1*3.8;

        for (size_t i = start; i < end; i++) {
                for (size_t j = 0; j < N; j++) {
                        if (j == i - 1 || j == i || j == i + 1)
                                continue;

                        if (protein_distance(self, i, j) < bead_diameter) {
                                dprintf("overlap of %g between atoms %u and %u\n",
                                        protein_distance(self, i, j), i, j);
                                return true;
                        }
                }
        }

        return false;
}



bool protein_do_end_move_first(struct protein *self, gsl_rng *rng)
{
        dprintf("moving first atom.\n");
        dprintf("before: atom(0) == "); dprint_vector(self->atom[0]);

        double R[3][3];
        make_random_rotation_matrix(R, rng);
        gsl_matrix_view RV = gsl_matrix_view_array((double *) R, 3, 3);
        rotate(false, &RV.matrix, self->atom[1], self->atom[0]);
        dprintf("after: atom(0) == "); dprint_vector(self->atom[0]);

        if (protein_is_overlapping(self, 0, 1)) {
                rotate(true, &RV.matrix, self->atom[1], self->atom[0]);
                dprintf("after undo: atom(0) == "); dprint_vector(self->atom[0]);
                return false;
        }

        return true;
}

bool protein_do_end_move_last(struct protein *self, gsl_rng *rng)
{
        dprintf("moving last atom.\n");
        dprintf("before: atom(%d) == ", self->num_atoms-1);
        dprint_vector(self->atom[self->num_atoms - 1]);

        const size_t N = self->num_atoms;

        double R[3][3];
        make_random_rotation_matrix(R, rng);
        gsl_matrix_view RV = gsl_matrix_view_array((double *) R, 3, 3);
        rotate(false, &RV.matrix, self->atom[N-2], self->atom[N-1]);

        dprintf("after: atom(%d) == ", N-1);
        dprint_vector(self->atom[N-1]);

        if (protein_is_overlapping(self, N-1, N)) {
                rotate(true, &RV.matrix,
                       self->atom[N - 2],
                       self->atom[N - 1]);
                dprintf("after undo: atom(%d) == ", N-1);
                dprint_vector(self->atom[N-1]);
                return false;
        }
        return true;
}



bool protein_do_shift_move(struct protein *self, gsl_rng *rng, size_t k)
{
        assert(k <= self->num_atoms - 3);

        bool status = true;

        dprintf("shifting move at atom %u\n", k);
        dprintf("before: atom(%u) == ", self->num_atoms-1); dprint_vector(self->atom[self->num_atoms-1]);
        dprintf("before: atom(%u) == ", k+1); dprint_vector(self->atom[k+1]);

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

        dprintf("after: atom(%u) == ", self->num_atoms-1); dprint_vector(self->atom[self->num_atoms-1]);
        dprintf("after: atom(%u) == ", k+1); dprint_vector(self->atom[k+1]);

        if (protein_is_overlapping(self, k+1, self->num_atoms)) {
                dprintf("undoing shift movement.\n");

                for (size_t i = k+2; i < self->num_atoms; i++)
                        gsl_vector_add(self->atom[i], t);

                dprintf("after undo: atom(%u) == ", self->num_atoms-1); dprint_vector(self->atom[self->num_atoms-1]);

                gsl_vector_memcpy(self->atom[k+1], bak);

                dprintf("after undo: atom(%u) == ", k+1); dprint_vector(self->atom[k+1]);

                status = false;
        }

        gsl_vector_free(bak);

        return status;
}



bool protein_do_spike_move(struct protein *self, gsl_rng *rng, size_t k)
{
        bool status = true;

        const double theta = 2*M_PI*gsl_rng_uniform_pos(rng);

        gsl_vector *bak = gsl_vector_alloc(3);
        gsl_vector_memcpy(bak, self->atom[k]);

        dprintf("before: atom(%u) == ", k); dprint_vector(self->atom[k]);

        /* v = p3-p1 */
        gsl_vector *v = gsl_vector_alloc(3);
        gsl_vector_memcpy(v, self->atom[k+1]);
        gsl_vector_sub(v, self->atom[k-1]);
        vector_normalize(v);

        /* a = p2-p1 */
        gsl_vector *a = gsl_vector_alloc(3);
        gsl_vector_memcpy(a, self->atom[k]);
        gsl_vector_sub(a, self->atom[k-1]);

        /* w = <a, v> v */
        gsl_vector *w = gsl_vector_alloc(3);
        gsl_vector_memcpy(w, v);
        double len;
        gsl_blas_ddot(a, v, &len);
        gsl_vector_scale(w, len);

        /* t = p1 + w */
        gsl_vector *t = gsl_vector_alloc(3);
        gsl_vector_memcpy(t, w);
        gsl_vector_add(t, self->atom[k-1]);

        /* q = atom[k] - t */
        gsl_vector *q = gsl_vector_alloc(3);
        gsl_vector_memcpy(q, self->atom[k]);
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
        assert(gsl_fcmp(gsl_blas_dnrm2(&u1.vector), 0.0, 1e-3) != 0);

        /* XXX This code could be replaced by BLAS' DROT. */
        double R[3][3] = {{cos(theta),  -sin(theta),    0.0},
                          {sin(theta),   cos(theta),    0.0},
                          {       0.0,          0.0,    1.0}};
        gsl_matrix_const_view RV = gsl_matrix_const_view_array((double *) R, 3, 3);

        gsl_matrix *A = gsl_matrix_alloc(3, 3);
        gsl_matrix *B = gsl_matrix_alloc(3, 3);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &RV.matrix, G, 0.0, A);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, G, A, 0.0, B);

        gsl_vector *z = gsl_vector_alloc(3);
        gsl_blas_dgemv(CblasNoTrans, 1.0, B, q, 0.0, z);

        gsl_vector_add(z, t);
        gsl_vector_memcpy(self->atom[k], z);

        dprintf("after: atom(%d) == ", k); dprint_vector(self->atom[k]);

        /* We are done if the conformation is correct. */
        if (protein_is_overlapping(self, k, k+1)) {
                gsl_vector_memcpy(self->atom[k], bak);
                dprintf("after undo: atom(%d) == ", k); dprint_vector(self->atom[k]);

                status = false;
        }

        gsl_vector_free(bak);
        gsl_vector_free(v); gsl_vector_free(a); gsl_vector_free(w);
        gsl_vector_free(t); gsl_vector_free(q); gsl_vector_free(z);
        gsl_matrix_free(G); gsl_matrix_free(A); gsl_matrix_free(B);

        return status;
}



bool protein_do_pivot_move(struct protein *self, gsl_rng *rng, size_t k)
{
        const double theta = 2*M_PI*gsl_rng_uniform_pos(rng);

        dprintf("pivoting move at atom %u.\n", k);
        dprintf("before: atom(%d) == ", k+1);
        dprint_vector(self->atom[k+1]);


        /* XXX This code could be replaced by BLAS' DROT. */
        double RR[3][3] = {{cos(theta), -sin(theta), 0.0},
                           {sin(theta),  cos(theta), 0.0},
                           {       0.0,         0.0, 1.0}};
        gsl_matrix_const_view RV = gsl_matrix_const_view_array((double *) RR, 3, 3);
        declare_stack_allocated_vector(y, 3);
        for (size_t i = k+1; i < self->num_atoms; i++) {
                gsl_vector_memcpy(y, self->atom[k]);
                gsl_vector_sub(self->atom[i], y);
                gsl_blas_dgemv(CblasNoTrans, 1.0, &RV.matrix, self->atom[i], 1.0, y);
                gsl_vector_memcpy(self->atom[i], y);
        }

        dprintf("after: atom(%d) == ", k+1);
        dprint_vector(self->atom[k+1]);

        if (protein_is_overlapping(self, k+1, self->num_atoms)) {
                dprintf("undoing invalid conformation.\n");
                for (size_t i = k+1; i < self->num_atoms; i++) {
                        gsl_vector_memcpy(y, self->atom[k]);
                        gsl_vector_sub(self->atom[i], y);
                        gsl_blas_dgemv(CblasTrans, 1.0,
                                       &RV.matrix, self->atom[i],
                                       1.0, y);
                        gsl_vector_memcpy(self->atom[i], y);
                }
                dprintf("after undo: atom(%d) == ", k+1);
                dprint_vector(self->atom[k+1]);

                return false;
        }

        return true;
}


void protein_scramble(struct protein *self, gsl_rng *rng)
{
        for (size_t r = 0; r < 10*self->num_atoms; r++)
                for (size_t k = 0; k < self->num_atoms; k++)
                        protein_do_natural_movement(self, rng, k);
}
