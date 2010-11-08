/**
 * @file protein.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "protein.h"
#include "utils.h"


const char gnuplot_command_line[] = GNUPLOT_EXECUTABLE " --persist";


struct protein *new_protein(size_t num_atoms, const double *positions)
{
        struct protein *m;

        assert(num_atoms != 0);
        assert(positions != NULL);

        if ((m = malloc(sizeof(struct protein))) == NULL)
                return NULL;

        m->num_atoms = num_atoms;

        m->positions = gsl_matrix_alloc(num_atoms, 3);
        if (m->positions == NULL) {
                free(m);
                return NULL;
        }

        gsl_matrix_const_view p = gsl_matrix_const_view_array(positions, num_atoms, 3);
        gsl_matrix_memcpy(m->positions, &p.matrix);

        m->native_contacts = gsl_matrix_calloc(num_atoms, num_atoms);
        protein_compute_contact_map(m, 10, m->native_contacts);

        return m;
}

#include "sample-proteins.c"

void delete_protein(struct protein *self)
{
        assert(self != NULL);
        assert(self->positions != NULL);
        gsl_matrix_free(self->positions);
        gsl_matrix_free(self->native_contacts);
        free(self);
}



gsl_vector_view protein_get_atom(const struct protein *self, size_t i)
{
        assert(self != NULL);
        assert(i < self->num_atoms);

        return gsl_matrix_row(self->positions, i);
}

static double signum(const struct protein *self, size_t i, size_t j)
{
        if (abs((int) (i - j)) == 3) {
                double uu[3], vv[3], ww[3];
                gsl_vector_view u = gsl_vector_view_array((double *) &uu, 3);
                gsl_vector_view v = gsl_vector_view_array((double *) &vv, 3);
                gsl_vector_view w = gsl_vector_view_array((double *) &ww, 3);

                gsl_vector_view v0 = gsl_matrix_row(self->positions, i);
                gsl_vector_view v1 = gsl_matrix_row(self->positions, i+1);
                gsl_vector_view v2 = gsl_matrix_row(self->positions, i+2);
                gsl_vector_view v3 = gsl_matrix_row(self->positions, i+3);

                gsl_vector_memcpy(&u.vector, &v1.vector);
                gsl_vector_memcpy(&v.vector, &v2.vector);
                gsl_vector_memcpy(&w.vector, &v3.vector);

                gsl_vector_sub(&u.vector, &v0.vector);
                gsl_vector_sub(&v.vector, &v1.vector);
                gsl_vector_sub(&w.vector, &v2.vector);

                double d = triple_scalar_product(&u.vector,
                                                 &v.vector,
                                                 &w.vector);
                return d/fabs(d);
        } else {
                return 1.0;
        }
}

double protein_distance(const struct protein *self, size_t i, size_t j)
{
        assert(self != NULL);

        gsl_vector_const_view v1 = gsl_matrix_const_row(self->positions, i);
        gsl_vector_const_view v2 = gsl_matrix_const_row(self->positions, j);

        double data[] = {0, 0, 0};
        gsl_vector_view uminusv = gsl_vector_view_array((double *) &data, 3);
        gsl_vector_memcpy(&uminusv.vector, &v1.vector);
        gsl_vector_sub(&uminusv.vector, &v2.vector);

        return signum(self, i, j)*gsl_blas_dnrm2(&uminusv.vector);
}

void protein_compute_contact_map(struct protein *self, double d_max, gsl_matrix *contact_map)
{
        assert(self != NULL);
        assert(d_max > 0);
        assert(contact_map != NULL);

        gsl_matrix_set_zero(contact_map);

        for (size_t i = 0; i < self->num_atoms; i++) {
                for (size_t j = i+1; j < self->num_atoms; j++) {
                        double d = protein_distance(self, i, j);
                        double b = (fabs(d) <= d_max) ? d : GSL_POSINF;

                        gsl_matrix_set(contact_map, i, j, b);
                        gsl_matrix_set(contact_map, j, i, b);
                }
        }
}



static int plot_structure(const gsl_matrix *positions)
{
        assert(positions != NULL);

        FILE *g = popen(gnuplot_command_line, "w");
        if (g == NULL)
                return -1;

        fprintf(g, "set terminal wxt\n"
                   "set view equal\n"
                   "set linetype 1 linecolor rgb 'dark-violet' linewidth 5\n"
                   "set title 'Alpha carbon structure'\n"
                   "unset xtics\n"
                   "unset ytics\n"
                   "unset ztics\n"
                   "unset border\n"
                   "splot '-' title '' with lines\n");

        print_matrix(g, positions);

        fprintf(g, "e\n");

        pclose(g);

        return 0;
}

static int plot_contact_map(const gsl_matrix *contact_map)
{
        assert(contact_map != NULL);

        FILE *g = popen(gnuplot_command_line, "w");
        if (g == NULL)
                return -1;

        fprintf(g, "set title 'Contact map'\n"
                   "set palette gray\n"
                   "unset colorbox\n"
                   "plot '-' matrix title '' with image\n");

        print_matrix(g, contact_map);

        fprintf(g, "e\n");

        pclose(g);

        return 0;
}

void protein_plot(const struct protein *self)
{
        plot_structure(self->positions);

        plot_contact_map(self->native_contacts);
}
