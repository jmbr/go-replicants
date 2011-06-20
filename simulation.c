#include "molecular-simulator.h"


const char U_file_template[] = "U--t-%02.05f--dmax-%02.05f--a-%02.05f.dat";
const char X_file_template[] = "X--t-%02.05f--dmax-%02.05f--a-%02.05f.xyz";


static int open_log_files(struct simulation *s);


struct simulation *new_simulation(const struct contact_map *native_map,
                                  double a, double temperature, gsl_rng *rng)
{
        if (native_map == NULL || a <= 0.0 || rng == NULL)
                return NULL;

        struct simulation *s = calloc(1, sizeof(struct simulation));
        if (s == NULL)
                return NULL;

        s->next_atom = 0;
        s->native_map = native_map;

        s->a = a;

        s->energy = GSL_POSINF;

        s->temperature = temperature;

        /* s->rng = rng; */
        s->rng = gsl_rng_alloc(gsl_rng_mt19937);
        struct timeval tv;
        gettimeofday(&tv, NULL);
        printf("Seeding the random number generator with %l\n", tv.tv_usec);
        gsl_rng_set(s->rng, tv.tv_usec);

        s->accepted = s->total = 0;

        if (open_log_files(s) == -1) {
                delete_simulation(s);
                return NULL;
        }

        return s;
}

int open_log_files(struct simulation *s)
{
        char name1[PATH_MAX], name2[PATH_MAX];

        sprintf(name1, X_file_template, s->temperature,
                contact_map_get_d_max(s->native_map), s->a);
        sprintf(name2, U_file_template, s->temperature,
                contact_map_get_d_max(s->native_map), s->a);

        s->X = fopen(name1, "a");
        s->U = fopen(name2, "a");

        return (s->X != NULL && s->U != NULL) ? 0 : -1;
}



void delete_simulation(struct simulation *self)
{
        assert(self);

        if (self->U != NULL)
                fclose(self->U);
        if (self->X != NULL)
                fclose(self->X);
        if (self->protein != NULL)
                delete_protein(self->protein);
        if (self->rng != NULL)
                gsl_rng_free(self->rng);

        free(self);
}


double simulation_get_acceptance_ratio(const struct simulation *self)
{
        assert(self != NULL);

        return (double) self->accepted/(double) self->total;
}


void simulation_first_iteration(struct simulation *self,
                                const struct protein *protein, double energy)
{
        assert(self != NULL);

        ++self->total;

        self->protein = protein_dup(protein);
        self->energy = energy;
}

static inline double
compute_potential_energy(const struct protein *x, const struct simulation *s)
{
        return potential(x, s->native_map, s->a);
}



void simulation_next_iteration(struct simulation *self)
{
        assert(self != NULL);

        ++self->total;

        struct protein *current = self->protein;
        struct protein *candidate = protein_dup(current);
        assert(candidate != NULL); /* XXX Add proper error checking here. */
        bool changed = protein_do_natural_movement(candidate, self->rng, self->next_atom);
        self->next_atom = (self->next_atom + 1) % self->protein->num_atoms;

        const double U1 = self->energy;
        const double U2 = changed ? compute_potential_energy(candidate, self) : U1;
        const double DU = U2 - U1;

        struct protein *chosen;

        if (DU <= 0.0)
                chosen = candidate;
        else {
                const double r = gsl_rng_uniform(self->rng);
                const double p = exp(-DU/self->temperature);
                chosen = r < p ? candidate : current;
        }

        if (chosen == candidate) {
                ++self->accepted;
                delete_protein(self->protein);
                self->protein = candidate;
                self->energy = U2;
        } else {
                delete_protein(candidate);
        }
}

void simulation_print_info(const struct simulation *self, FILE *stream)
{
        assert(self != NULL);

        fprintf(stream, "simulation (T = %02.2f, a = %2.1f, d_max = %2.1f)\n",
                self->temperature, self->a,
                contact_map_get_d_max(self->native_map));
}
