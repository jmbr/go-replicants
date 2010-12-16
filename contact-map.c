#include "molecular-simulator.h"


struct contact_map {
        size_t num_atoms;
        size_t num_contacts;
        double d_max;
        double *distance;
};


static void contact_map_compute(struct contact_map *self, const struct protein *protein);
static void print_contact_map(FILE *stream, const struct contact_map *self);


struct contact_map *new_contact_map(const struct protein *protein, double d_max)
{
        if (protein == NULL || d_max <= 0.0)
                return NULL;

        struct contact_map *c;

        if ((c = malloc(sizeof(struct contact_map))) == NULL)
                return NULL;

        const size_t N = protein->num_atoms;

        c->num_atoms = N;
        c->num_contacts = 0;
        c->d_max = d_max;
        c->distance = calloc(N*N, sizeof(double));
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

        if (self->distance)
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
                for (size_t j = i+2; j < N; j++) {
                        double d = protein_distance(protein, i, j);

                        switch (j - i) {
                        case 2:
                                self->distance[N*i + j] = d;
                                /* printf("%u, %u: %f\n", i, j, d); */
                                ++self->num_contacts;
                                break;
                        case 3:
                                self->distance[N*i + j] = protein_signum(protein, i, j)*d;
                                /* printf("%u, %u: %f\n", i, j, d); */
                                ++self->num_contacts;
                                break;
                        default:
                                if (d <= self->d_max) {
                                        self->distance[N*i + j] = d;
                                        /* printf("%u, %u: %f\n", i, j, d); */
                                        ++self->num_contacts;
                                }
                                break;
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

        const size_t N = self->num_atoms;

        if (i == j)
                return 0.0;
        else if (i < j)
                return self->distance[N*i + j];
        else
                return self->distance[N*j + i];
}



void contact_map_plot(const struct contact_map *self, FILE *gnuplot,
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
                         "set palette gray\n"
                         "unset colorbox\n"
                         "set tics out nooffset\n"
                         "set xtics 0,5,%u\n"
                         "set ytics 0,5,%u\n"
                         "plot '-' matrix notitle with image\n",
                self->num_atoms, self->num_atoms);
        print_contact_map(gnuplot, self);
        fprintf(gnuplot, "e\ne\n");

        fflush(gnuplot);
}

void print_contact_map(FILE *stream, const struct contact_map *self)
{
        const size_t N = self->num_atoms;

        for (size_t i = 0; i < N; i++) {
                for (size_t j = 0; j < N; j++) {
                        double d = contact_map_get_distance(self, i, j);
                        if (abs((int) i - (int) j) < 2)
                                fputs("0 ", stream);
                        else
                                fputs(d != 0.0 ? "0 " : "1 ", stream);
                }
                fputs("\n", stream);
        }
}

int contact_map_diff(const struct contact_map *m1, const struct contact_map *m2)
{
        if (m1->num_atoms != m2->num_atoms)
                return -1;

        const size_t N = m1->num_atoms;

        for (size_t i = 0; i < N; i++) {
                for (size_t j = i+2; j < N; j++) {
                        double d1 = contact_map_get_distance(m1, i, j);
                        double d2 = contact_map_get_distance(m2, i, j);

                        if (d1 != 0.0 && d2 == 0.0)
                                printf("contact between %u and %u is not present in the second map\n", i, j);
                        if (d1 == 0.0 && d2 != 0.0)
                                printf("contact between %u and %u is not present in the first map\n", i, j);
                }
        }

        return 0;
}
