#include "molecular-simulator.h"


const size_t simulation_save_step = 10000;

const char gnuplot_command_line[] = GNUPLOT_EXECUTABLE " -persist";

const size_t bufsize = 1024;
const char configurations_file_template[] = "configurations--t-%02.02f--dmax-%02.02f--a-%02.02f.dat";
const char energies_file_template[] = "energies--t-%02.02f--dmax-%02.02f--a-%02.02f.dat";


struct simulation *new_simulation(const struct protein *p,
                                  gsl_rng *rng,
                                  double temperature,
                                  const struct simulation_options *opts)
{
        if (p == NULL || rng == NULL || opts == NULL
            || opts->d_max <= 0.0 || opts->a <= 0.0)
        {
                return NULL;
        }

        struct simulation *s = calloc(1, sizeof(struct simulation));
        if (s == NULL)
                return NULL;

        s->protein = protein_dup(p);

        /* Initialize native contact map. */
        s->d_max = opts->d_max;
        if ((s->native_map = new_contact_map(s->protein, s->d_max)) == NULL)
                delete_simulation(s);

        s->a = opts->a;

        s->T = temperature;

        /* Compute the potential energy for the protein. */
        s->orig_energy = s->energy = potential(s->protein, s->native_map, s->a);

        s->rng = rng;

        if ((s->gnuplot = popen(gnuplot_command_line, "w")) == NULL)
                delete_simulation(s);

        s->accepted = 0;
        s->total = 0;

        char name1[bufsize], name2[bufsize];
        sprintf(name1, configurations_file_template, s->T, s->d_max, s->a);
        sprintf(name2, energies_file_template, s->T, s->d_max, s->a);
        s->configurations = fopen(name1, "w");
        s->energies = fopen(name2, "w");
        if (s->configurations == NULL || s->energies == NULL) {
                delete_simulation(s);
                return NULL;
        }

        s->next_atom = 1;

        return s;
}

void delete_simulation(struct simulation *self)
{
        assert(self);

        if (self->gnuplot != NULL)
                pclose(self->gnuplot);
        if (self->native_map != NULL)
                delete_contact_map(self->native_map);
        if (self->configurations != NULL)
                fclose(self->configurations);
        if (self->energies != NULL)
                fclose(self->energies);
        if (self->protein != NULL)
                delete_protein(self->protein);

        free(self);
}



double simulation_get_acceptance_ratio(const struct simulation *self)
{
        assert(self != NULL);

        return (double) self->accepted/(double) self->total;
}



static inline double
compute_potential_energy(const struct protein *x, const struct simulation *s)
{
        return potential(x, s->native_map, s->a);
}

static void simulation_save_state(const struct simulation *self)
{
        protein_print_atoms(self->protein, self->configurations);
        fflush(self->configurations);

        fprintf(self->energies, "%e\n", self->energy);
        fflush(self->energies);
}

void simulation_first_iteration(struct simulation *self)
{
        assert(self != NULL);

        ++self->total;

        protein_scramble(self->protein, self->rng);

        self->energy = compute_potential_energy(self->protein, self);

        simulation_save_state(self);
}

void simulation_next_iteration(struct simulation *self)
{
        assert(self != NULL);

        ++self->total;

        if (self->total % simulation_save_step == 0)
                simulation_save_state(self);

        struct protein *current = self->protein;
        struct protein *candidate = protein_dup(current);
        assert(candidate != NULL); /* XXX Check this properly. */

        protein_do_natural_movement(candidate, self->rng, self->next_atom);
        self->next_atom = (self->next_atom + 1) % self->protein->num_atoms;

        const double U1 = compute_potential_energy(current, self);
        const double U2 = compute_potential_energy(candidate, self);
        const double DU = U2 - U1;

        struct protein *chosen;

        if (DU <= 0.0)
                chosen = candidate;
        else {
                const double r = gsl_rng_uniform(self->rng);
                const double p = -1.0/self->T*DU;
                chosen = r < gsl_min(1.0, exp(p)) ? candidate : current;
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

bool simulation_has_converged(const struct simulation *self)
{
        assert(self != NULL);

        return gsl_fcmp(self->orig_energy, self->energy, 1e-3) == 0;
}



void simulation_print_info(const struct simulation *self, FILE *stream)
{
        assert(self != NULL);

        fprintf(stream, "simulation (T = %02.2f, a = %2.1f, d_max = %2.1f)\n",
                self->T, self->a, self->d_max);
}
