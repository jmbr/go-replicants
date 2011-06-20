#include "molecular-simulator.h"


const size_t save_energy_step = 1000;
const size_t save_conformation_step = 5000;

static bool options_are_invalid(const struct simulation_options *options);
static void replicas_exchange(struct replicas *self, size_t k);
static void save_energy(const struct replicas *self);
static void save_conformation(const struct replicas *self);

/* XXX Replicas should be responsible for allocating and freeing the
 * protein structure.  */
struct replicas *new_replicas(struct protein *protein,
                              const struct simulation_options *options)
{
        if (protein == NULL || options_are_invalid(options))
                return NULL;

        const size_t total_size = sizeof(struct replicas)
                + options->num_replicas*sizeof(struct simulation *);
        
        struct replicas *r = calloc(1, total_size);
        if (r == NULL)
                return NULL;

        r->rng = options->rng;
        r->protein = protein;
        r->native_map = new_contact_map(protein, options->d_max);
        if (r->native_map == NULL) {
                free(r);
                return NULL;
        }
        r->a = options->a;
        r->num_replicas = options->num_replicas;
        r->exchanges = calloc(r->num_replicas, sizeof(size_t));
        r->total = calloc(r->num_replicas, sizeof(size_t));
        if (r->exchanges == NULL || r->total == NULL) {
                delete_replicas(r);
                return NULL;
        }
        if ((r->log = fopen("replicas.log", "a")) == NULL)
                delete_replicas(r);
        for (size_t k = 0; k < r->num_replicas; k++) {
                r->replica[k] = new_simulation(r->native_map, r->a,
                                               options->temperatures[k], r->rng);
                if (r->replica[k] == NULL) {
                        delete_replicas(r);
                        return NULL;
                }
        }

        return r;
}

bool options_are_invalid(const struct simulation_options *options)
{
        return (options == NULL
                || options->rng == NULL || options->num_replicas == 0
                || options->d_max <= 0.0 || options->a <= 0.0);
}



void delete_replicas(struct replicas *self)
{
        for (size_t k = 0; k < self->num_replicas; k++)
                if (self->replica[k] != NULL)
                        delete_simulation(self->replica[k]);
        if (self->exchanges != NULL)
                free(self->exchanges);
        if (self->total != NULL)
                free(self->total);
        if (self->native_map != NULL)
                delete_contact_map(self->native_map);
        /* XXX: If new_replica fails, protein will be deallocated. */
        if (self->protein != NULL)
                delete_protein(self->protein);
        if (self->log != NULL)
                fclose(self->log);
        free(self);
}

void save_energy(const struct replicas *self)
{
        for (size_t k = 0; k < self->num_replicas; k++) {
                const struct simulation *s = self->replica[k];
                fprintf(s->U, "%f\n", s->energy);
                fflush(s->U);
        }
}

void save_conformation(const struct replicas *self)
{
        for (size_t k = 0; k < self->num_replicas; k++) {
                const struct simulation *s = self->replica[k];
                protein_write_xyz(s->protein, s->X);
        }
}

void replicas_first_iteration(struct replicas *self)
{
        protein_scramble(self->protein, self->rng);
        double energy = potential(self->protein, self->native_map, self->a);

        size_t k;
#pragma omp parallel for private(k)
        for (k = 0; k < self->num_replicas; k++)
                simulation_first_iteration(self->replica[k],
                                           self->protein, energy);

        save_energy(self);
        save_conformation(self);
}

void replicas_resume(struct replicas *self, const struct protein *conf[])
{
        /* Initialize every replica with the protein structure and energy. */
        size_t k;
#pragma omp parallel for private(k)
        for (k = 0; k < self->num_replicas; k++) {
                const struct protein *p = conf[k];
                const double U = potential(p, self->native_map, self->a);
                simulation_first_iteration(self->replica[k], p, U);
        }
}



void replicas_thermalize(struct replicas *self, size_t num_iters)
{
        const size_t iters_per_cycle = self->protein->num_atoms;

        fprintf(self->log, "performing %u thermalization steps.\n", num_iters);
        
        for (size_t s = 0; s < num_iters; s++) {
                size_t k;
#pragma omp parallel for private(k)
                for (k = 0; k < self->num_replicas; k++)
                        for (size_t c = 0; c < iters_per_cycle; c++)
                                simulation_next_iteration(self->replica[k]);

                if (s > 0 && s % save_energy_step == 0)
                        save_energy(self);
                if (s > 0 && s % save_conformation_step == 0)
                        save_conformation(self);
        }

        fprintf(self->log, "done with the thermalization steps.\n");
}

void replicas_next_iteration(struct replicas *self)
{
        size_t k;

        const size_t num_iters = 5000;

        if (self->num_replicas > 1) {
                fprintf(self->log, "attempting to exchange replicas.\n");
                for (k = gsl_rng_uniform_int(self->rng, 2);
                     k <= self->num_replicas - 2;
                     k += 2)
                {
                        replicas_exchange(self, k);
                }
                fprintf(self->log, "done with replica exchange.\n");
        }

        for (size_t s = 0; s < num_iters; s++) {
#pragma omp parallel for private(k)
                for (k = 0; k < self->num_replicas; k++)
                        for (size_t c = 0; c < self->protein->num_atoms; c++)
                                simulation_next_iteration(self->replica[k]);

                if (s % save_energy_step == 0)
                        save_energy(self);
                if (s % save_conformation_step == 0)
                        save_conformation(self);

                fflush(self->log);
        }
}



void replicas_exchange(struct replicas *self, size_t k)
{
        assert(self->num_replicas >= 2);

        struct simulation *s1 = self->replica[k];
        struct simulation *s2 = self->replica[k+1];
        const double U1 = s1->energy;
        const double U2 = s2->energy;
        const double T1 = s1->temperature;
        const double T2 = s2->temperature;
        const double DU = U2 - U1;
        const double DB = 1.0/T2 - 1.0/T1;
        const double p = exp(DB*DU);
        const double r = gsl_rng_uniform(self->rng);
        if (r < p) {
                fprintf(self->log, "swapping replicas %u and %u.\n", k, k+1);
                struct protein *x = s1->protein;
                s1->protein = s2->protein;
                s1->energy = U2;
                s2->protein = x;
                s2->energy = U1;
                ++self->exchanges[k];
        }

        ++self->total[k];
}



size_t replicas_total_exchanges(const struct replicas *self)
{
        size_t total = 0;
        for (size_t r = 0; r < self->num_replicas - 1; r++)
                total += self->exchanges[r];
        return total;
}

void replicas_get_exchange_ratios(const struct replicas *self,
                                  double ratios[])
{
        for (size_t r = 0; r < self->num_replicas - 1; r++)
                ratios[r] = (double) self->exchanges[r]/self->total[r];
}

void replicas_print_info(const struct replicas *self, FILE *stream)
{
        double ratios[self->num_replicas - 1];

        replicas_get_exchange_ratios(self, ratios);

        for (size_t r = 0; r < self->num_replicas - 1; r++)
                fprintf(stream, "ratio of exchanges between temperatures "
                        "%2.2f and %2.2f: %2.2f\n",
                        self->replica[r]->temperature,
                        self->replica[r+1]->temperature, ratios[r]);
        fflush(stream);
}
