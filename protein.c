#include "molecular-simulator.h"


static struct protein *read_xyz_file(FILE *f);
static int write_xyz_file(const struct protein *self, FILE *stream);

static void protein_do_shift_move(struct protein *self, gsl_rng *rng, size_t k,
                                  bool undo);
static void protein_do_spike_move(struct protein *self, gsl_rng *rng, size_t k,
                                  bool undo);
static void protein_do_pivot_move(struct protein *self, gsl_rng *rng, size_t k,
                                  bool undo);
static void protein_do_end_move_first(struct protein *self, gsl_rng *rng, bool undo);
static void protein_do_end_move_last(struct protein *self, gsl_rng *rng, bool undo);


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

struct protein *protein_dup(const struct protein *self)
{
        if (self == NULL)
                return NULL;

        struct protein *p = malloc(sizeof(struct protein)
                                   + self->num_atoms*sizeof(gsl_vector *));

        assert(p != NULL);

        /* XXX This can be sped up by allocating a large gsl_block
         * beforehand.
         */
        p->num_atoms = self->num_atoms;
        for (size_t i = 0; i < p->num_atoms; i++) {
                p->atom[i] = gsl_vector_alloc(3);
                gsl_vector_memcpy(p->atom[i], self->atom[i]);
        }

        return p;
}



struct protein *protein_read_xyz_file(const char *name)
{
        FILE *f;

        if ((f = fopen(name, "r")) == NULL)
                return NULL;

        struct protein *p = read_xyz_file(f);

        fclose(f);

        return p;
}

struct protein *read_xyz_file(FILE *f)
{
        size_t num_atoms;

        if (fscanf(f, "%u\n", &num_atoms) != 1)
                return NULL;

        if (fscanf(f, "%*s") == EOF)
                return NULL;

        double *atom_tab = calloc(num_atoms, 3*sizeof(double));

        for (size_t k = 0; k < num_atoms; k++) {
                int status;
                status = fscanf(f, "%*s %lg %lg %lg",
                                atom_tab + 3*k + 0,
                                atom_tab + 3*k + 1,
                                atom_tab + 3*k + 2);
                if (status != 3) {
                        free(atom_tab);
                        return NULL;
                }
        }

        return new_protein(num_atoms, atom_tab);
}



int protein_write_xyz_file(const struct protein *self, const char *name)
{
        FILE *f;

        if ((f = fopen(name, "w")) == NULL)
                return -1;

        int status = write_xyz_file(self, f);

        fclose(f);

        return status;
}

int write_xyz_file(const struct protein *self, FILE *stream)
{
        int status, n = 0;

        fprintf(stream, "%u\n", self->num_atoms);
        fprintf(stream, "Protein\n");
        for (size_t i = 0; i < self->num_atoms; i++) {
                status = fprintf(stream, "CA %g %g %g\n",
                                 gsl_vector_get(self->atom[i], 0),
                                 gsl_vector_get(self->atom[i], 1),
                                 gsl_vector_get(self->atom[i], 2));
                if (status < 0)
                        return -1;
                n += status;
        }

        return n;
}



void protein_plot(const struct protein *self, FILE *gnuplot,
                  const char *title_format, ...)
{
        assert(self != NULL);
        assert(gnuplot != NULL);

        fprintf(gnuplot, "set terminal wxt 0 noraise\n");

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
                status = fprintf(stream, "%f %f %f\n",
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



double protein_signum(const struct protein *self, size_t i, size_t j)
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

        return gsl_blas_dnrm2(&v.vector);
}



void protein_do_movement(struct protein *self, gsl_rng *rng,
                         enum protein_movements m, size_t k, bool undo)
{
        assert(self != NULL);
        assert(rng != NULL);

        switch (m) {
        case PROTEIN_SPIKE_MOVE:
                protein_do_spike_move(self, rng, k, undo);
                break;
        case PROTEIN_SHIFT_MOVE:
                protein_do_shift_move(self, rng, k, undo);
                break;
        case PROTEIN_PIVOT_MOVE:
                protein_do_pivot_move(self, rng, k, undo);
                break;
        case PROTEIN_END_MOVE_FIRST:
                protein_do_end_move_first(self, rng, undo);
                break;
        case PROTEIN_END_MOVE_LAST:
                protein_do_end_move_last(self, rng, undo);
                break;
        default:
                assert(0);
                break;
        }
}

void protein_do_natural_movement(struct protein *self, gsl_rng *rng, size_t k)
{
        dprintf("thread #%d is changing atom %u\n", omp_get_thread_num(), k);

        if (k == 1)
                protein_do_end_move_first(self, rng, true);
        else if (1 < k && k < self->num_atoms-2)
                protein_do_movement(self, rng, (enum protein_movements) gsl_rng_uniform_int(rng, 3), k, true);
        else
                protein_do_end_move_last(self, rng, true);
}



void protein_do_end_move_first(struct protein *self, gsl_rng *rng, bool undo)
{
        dprintf("moving first atom.\n");
        dprintf("before: atom(0) == "); dprint_vector(self->atom[0]);

        double R[3][3];
        make_random_rotation_matrix(R, rng);
        gsl_matrix_view RV = gsl_matrix_view_array((double *) R, 3, 3);
        rotate(false, &RV.matrix, self->atom[1], self->atom[0]);
        dprintf("after: atom(0) == "); dprint_vector(self->atom[0]);

        if (undo && protein_is_overlapping(self)) {
                rotate(true, &RV.matrix, self->atom[1], self->atom[0]);
                dprintf("after undo: atom(0) == "); dprint_vector(self->atom[0]);
        }
}

void protein_do_end_move_last(struct protein *self, gsl_rng *rng, bool undo)
{
        dprintf("moving last atom.\n");
        dprintf("before: atom(%d) == ", self->num_atoms-1);
        dprint_vector(self->atom[self->num_atoms - 1]);

        double R[3][3];
        make_random_rotation_matrix(R, rng);
        gsl_matrix_view RV = gsl_matrix_view_array((double *) R, 3, 3);
        rotate(false, &RV.matrix, self->atom[self->num_atoms - 2],
                           self->atom[self->num_atoms - 1]);

        dprintf("after: atom(%d) == ", self->num_atoms-1);
        dprint_vector(self->atom[self->num_atoms - 1]);

        if (undo && protein_is_overlapping(self)) {
                rotate(true, &RV.matrix,
                       self->atom[self->num_atoms - 2],
                       self->atom[self->num_atoms - 1]);
                dprintf("after undo: atom(%d) == ", self->num_atoms-1);
                dprint_vector(self->atom[self->num_atoms-1]);
        }
}



void protein_do_shift_move(struct protein *self, gsl_rng *rng, size_t k, bool undo)
{
        assert(k <= self->num_atoms - 3);

        /* const size_t k = gsl_rng_uniform_int(rng, self->num_atoms - 3); */

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

        if (undo && protein_is_overlapping(self)) {
                dprintf("undoing shift movement.\n");

                for (size_t i = k+2; i < self->num_atoms; i++)
                        gsl_vector_add(self->atom[i], t);

                dprintf("after undo: atom(%u) == ", self->num_atoms-1); dprint_vector(self->atom[self->num_atoms-1]);

                gsl_vector_memcpy(self->atom[k+1], bak);

                dprintf("after undo: atom(%u) == ", k+1); dprint_vector(self->atom[k+1]);
        }

        gsl_vector_free(bak);
}



void protein_do_spike_move(struct protein *self, gsl_rng *rng, size_t k, bool undo)
{
        assert(k <= self->num_atoms - 3);

        const double theta = 2*M_PI*gsl_rng_uniform_pos(rng);

        gsl_vector *bak = gsl_vector_alloc(3);
        gsl_vector_memcpy(bak, self->atom[k+1]);

        dprintf("before: atom(%u) == ", k+1); dprint_vector(self->atom[k+1]);

        /* v = p3-p1 */
        gsl_vector *v = gsl_vector_alloc(3);
        gsl_vector_memcpy(v, self->atom[k+2]);
        gsl_vector_sub(v, self->atom[k]);
        vector_normalize(v);

        /* a = p2-p1 */
        gsl_vector *a = gsl_vector_alloc(3);
        gsl_vector_memcpy(a, self->atom[k+1]);
        gsl_vector_sub(a, self->atom[k]);

        /* w = <a, v> v */
        gsl_vector *w = gsl_vector_alloc(3);
        gsl_vector_memcpy(w, v);
        double len;
        gsl_blas_ddot(a, v, &len);
        gsl_vector_scale(w, len);

        /* t = p1 + w */
        gsl_vector *t = gsl_vector_alloc(3);
        gsl_vector_memcpy(t, w);
        gsl_vector_add(t, self->atom[k]);

        /* q = atom[k+1] - t */
        gsl_vector *q = gsl_vector_alloc(3);
        gsl_vector_memcpy(q, self->atom[k+1]);
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
        gsl_vector_memcpy(self->atom[k+1], z);

        dprintf("after: atom(%d) == ", k+1); dprint_vector(self->atom[k+1]);

        /* We are done if the conformation is correct. */
        if (undo && protein_is_overlapping(self)) {
                gsl_vector_memcpy(self->atom[k+1], bak);
                dprintf("after undo: atom(%d) == ", k+1); dprint_vector(self->atom[k+1]);
        }

        gsl_vector_free(bak);
        gsl_vector_free(v); gsl_vector_free(a); gsl_vector_free(w);
        gsl_vector_free(t); gsl_vector_free(q); gsl_vector_free(z);
        gsl_matrix_free(G); gsl_matrix_free(A); gsl_matrix_free(B);
}



void protein_do_pivot_move(struct protein *self, gsl_rng *rng, size_t k, bool undo)
{
        const double theta = 2*M_PI*gsl_rng_uniform_pos(rng);

        dprintf("pivoting move at atom %u.\n", k);
        dprintf("before: atom(%d) == ", k+1);
        dprint_vector(self->atom[k+1]);

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
        dprint_vector(self->atom[k+1]);

        if (undo && protein_is_overlapping(self)) {
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
        }
}



void protein_scramble(struct protein *self, gsl_rng *rng)
{
        for (size_t r = 0; r < 10*self->num_atoms; r++)
                for (size_t k = 0; k < self->num_atoms; k++)
                        protein_do_natural_movement(self, rng, k);
}

