/**
 * @file contact_map.c
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <assert.h>

#include <gsl/gsl_math.h>

#include "protein.h"
#include "contact-map.h"


struct contact_map {
        size_t num_atoms;
        size_t num_contacts;
        double d_max;
        double *distance;
};


static void contact_map_compute(struct contact_map *self, const struct protein *protein);
static int print_contact_map(FILE *stream, const struct contact_map *self);
static inline size_t idx(size_t N, size_t i, size_t j)
{
        if (i == 0) /* XXX Is this the most compact expression possible? */
                return j - i - 1;
        else
                return j - i - 1 + i*(N-1) - (i-1)*i/2;
}


        
struct contact_map *new_contact_map(const struct protein *protein, double d_max)
{
        if (protein == NULL || d_max <= 0.0)
                return NULL;

        struct contact_map *c;

        if ((c = malloc(sizeof(struct contact_map))) == NULL)
                return NULL;

        c->num_atoms = protein->num_atoms;
        c->num_contacts = 0;
        c->d_max = d_max;
        c->distance = calloc(c->num_atoms*(c->num_atoms - 1)/2, sizeof(double));
        if (c->distance == NULL) {
                free(c);
                return NULL;
        }

        contact_map_compute(c, protein);

        return c;
}

void delete_contact_map(struct contact_map *self)
{
        assert(self != NULL);
        assert(self->distance != NULL);

        free(self->distance);
        free(self);
}

void contact_map_compute(struct contact_map *self,
                         const struct protein *protein)
{
        assert(self != NULL);
        assert(protein != NULL);

        const size_t N = self->num_atoms;

        for (size_t i = 0; i < N; i++) {
                for (size_t j = i+1; j < N; j++) {
                        double d = protein_distance(protein, i, j);

                        if (d <= self->d_max) {
                                self->distance[idx(N, i, j)] = protein_signum(protein, i, j)*d;
                                ++self->num_contacts;
                        } else {
                                self->distance[idx(N, i, j)] = GSL_POSINF;
                        }
                }
        }
}



size_t contact_map_get_num_contacts(const struct contact_map *self)
{
        assert(self != NULL);

        return self->num_contacts;
}

double contact_map_get_d_max(const struct contact_map *self)
{
        assert(self != NULL);

        return self->d_max;
}

double contact_map_get_distance(const struct contact_map *self, 
                                size_t i, size_t j)
{
        assert(self != NULL);
        assert(self->num_atoms > 0);
        assert(i < self->num_atoms);
        assert(j < self->num_atoms);

        if (i == j)
                return 0.0;
        else if (i < j)
                return self->distance[idx(self->num_atoms, i, j)];
        else
                return self->distance[idx(self->num_atoms, j, i)];
}



int contact_map_plot(const struct contact_map *self, FILE *gnuplot)
{
        assert(self != NULL);

        if (gnuplot == NULL)
                return -1;

        int status = fprintf(gnuplot,
                             "set title 'Contact map ($d_max$ = %2.3g)'\n"
                             "set palette gray\n"
                             "unset colorbox\n"
                             "plot '-' matrix title '' with image\n",
                             self->d_max);
        if (status < 0)
                return -1;

        if ((status = print_contact_map(gnuplot, self)) < 0)
                return -1;

        if (fprintf(gnuplot, "e\ne\n") < 0)
                return -1;

        if (fflush(gnuplot) != 0)
                return -1;

        return 0;
}

int print_contact_map(FILE *stream, const struct contact_map *self)
{
        int status, n = 0;
        const size_t N = self->num_atoms;

        for (size_t i = 0; i < N; i++) {
                for (size_t j = 0; j < N; j++) {
                        double d = contact_map_get_distance(self, i, j);
                        status = fputs(d <= self->d_max ? "0 " : "1 ", stream);
                        if (status < 0)
                                return -1;
                        n += status;
                }

                if (fprintf(stream, "\n") < 0)
                        return -1;
                ++n;
        }

        return n;
}
