#include "molecular-simulator.h"


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

        struct protein *p = protein_read_xyz(f);

        fclose(f);

        return p;
}

struct protein *protein_read_xyz(FILE *stream)
{
        size_t num_atoms;

        if (fscanf(stream, "%u\n", &num_atoms) != 1)
                return NULL;

        if (fscanf(stream, "%*s") == EOF)
                return NULL;

        double *atom_tab = calloc(num_atoms, 3*sizeof(double));

        for (size_t k = 0; k < num_atoms; k++) {
                int status;
                status = fscanf(stream, "%*s %lg %lg %lg",
                                atom_tab + 3*k + 0,
                                atom_tab + 3*k + 1,
                                atom_tab + 3*k + 2);
                if (status != 3) {
                        free(atom_tab);
                        return NULL;
                }
        }

        struct protein *p = new_protein(num_atoms, atom_tab);
        free(atom_tab);

        return p;
}

struct protein *protein_read_latest_xyz(FILE *stream)
{
        if (stream == NULL)
                return NULL;

        struct protein *p, *q = NULL;

        for (p = protein_read_xyz(stream);
             p != NULL;
             p = protein_read_xyz(stream))
        {
                if (q)
                        delete_protein(q);
                q = p;
        }

        return q;
}



int protein_write_xyz_file(const struct protein *self, const char *name)
{
        FILE *f;

        if ((f = fopen(name, "w")) == NULL)
                return -1;

        int status = protein_write_xyz(self, f);

        fclose(f);

        return status;
}

int protein_write_xyz(const struct protein *self, FILE *stream)
{
        int status, n = 0;

        fprintf(stream, "%u\n"
                        "Protein\n", self->num_atoms);
        for (size_t i = 0; i < self->num_atoms; i++) {
                status = fprintf(stream, "CA %g %g %g\n",
                                 gsl_vector_get(self->atom[i], 0),
                                 gsl_vector_get(self->atom[i], 1),
                                 gsl_vector_get(self->atom[i], 2));
                if (status < 0)
                        return -1;
                n += status;
        }

        fflush(stream);

        return n;
}



void protein_plot(const struct protein *self, FILE *gnuplot, bool draw_labels,
                  const char *title_format, ...)
{
        assert(self != NULL);
        assert(gnuplot != NULL);

        va_list ap;
        fprintf(gnuplot, "set title '");
        va_start(ap, title_format);
        vfprintf(gnuplot, title_format, ap);
        va_end(ap);
        fprintf(gnuplot, "'\n");

        fprintf(gnuplot, "set terminal wxt noraise\n"
                         "set view equal xyz\n"
                         "set linetype 1 linecolor palette z linewidth 5\n"
                         "unset tics\n"
                         "unset border\n"
                         "unset colorbox\n"
                         "splot '-' notitle with lines\n");

        protein_print_atoms(self, gnuplot);

        fprintf(gnuplot, "e\n");

        if (draw_labels)
                for (size_t i = 0; i < self->num_atoms; i++)
                        fprintf(gnuplot, "set label '%u' at %f, %f, %f\n", i,
                                gsl_vector_get(self->atom[i], 0),
                                gsl_vector_get(self->atom[i], 1),
                                gsl_vector_get(self->atom[i], 2));

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



double protein_signum(const struct protein *self, size_t i, size_t j)
{
        if (abs((int) (i - j)) == 3) { /* XXX We don't need this check if we make sure j > i */
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

                /* XXX Should we check that this is not zero? */
                return signbit(triple_scalar_product(&u.vector, &v.vector, &w.vector)) != 0 ? -1.0 : 1.0;
        } else {
                return 1.0;
        }
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
