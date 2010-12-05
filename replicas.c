#include "molecular-simulator.h"


static void replicas_exchange(struct replicas *self, size_t k);


struct replicas *new_replicas(struct protein *protein,
                              const struct simulation_options *options)
{
        if (protein == NULL
            || options == NULL
            || options->rng == NULL
            || options->num_replicas == 0
            || options->d_max <= 0.0 || options->a <= 0.0)
        {
                return NULL;
        }
        
        struct replicas *r;
        r = calloc(1, sizeof(struct replicas)
                      + options->num_replicas*sizeof(struct simulation *));
        if (r == NULL)
                return NULL;

        r->protein = protein;

        r->rng = options->rng;

        /* Initialize native contact map. */
        r->native_map = new_contact_map(protein, options->d_max);
        if (r->native_map == NULL) {
                free(r);
                return NULL;
        }

        r->a = options->a;

        r->orig_energy = potential(r->protein, r->native_map, r->a);

        r->exchanges = calloc(options->num_replicas, sizeof(size_t));
        r->total = calloc(options->num_replicas, sizeof(size_t));
        if (r->exchanges == NULL || r->total == NULL) {
                delete_replicas(r);
                return NULL;
        }

        r->num_replicas = options->num_replicas;

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

void delete_replicas(struct replicas *self)
{
        for (size_t k = 0; k < self->num_replicas; k++)
                if (self->replica[k] != NULL)
                        delete_simulation(self->replica[k]);
        if (self->exchanges)
                free(self->exchanges);
        if (self->total)
                free(self->total);
        if (self->native_map)
                delete_contact_map(self->native_map);
        if (self->protein)
                delete_protein(self->protein);
        free(self);
}



void replicas_first_iteration(struct replicas *self)
{
        protein_scramble(self->protein, self->rng);
        double energy = potential(self->protein, self->native_map, self->a);

        protein_write_xyz_file(self->protein, "initial-conformation.xyz");

        size_t k;
        #pragma omp parallel for private(k)
        for (k = 0; k < self->num_replicas; k++)
                simulation_first_iteration(self->replica[k],
                                           self->protein, energy);
}

void replicas_resume(struct replicas *self, const struct protein *config[])
{
        size_t k;
        #pragma omp parallel for private(k)
        for (k = 0; k < self->num_replicas; k++) {
                const struct protein *p = config[k];
                const double U = potential(p, self->native_map, self->a);
                simulation_first_iteration(self->replica[k], p, U);
        }
}

void replicas_thermalize(struct replicas *self, size_t num_iters)
{
        const size_t iters_per_cycle = self->protein->num_atoms;

        dprintf("performing %u thermalization steps.\n", num_iters);

        replicas_first_iteration(self);

        for (size_t s = 0; s < (num_iters-1)*iters_per_cycle; s++) {
                size_t k;
                #pragma omp parallel for private(k)
                for (k = 0; k < self->num_replicas; k++)
                        simulation_next_iteration(self->replica[k]);
        }

        dprintf("done with the thermalization steps.\n");
}

void replicas_next_iteration(struct replicas *self)
{
        size_t k;

        const size_t num_iters = 5000;
        const size_t iters_per_cycle = self->protein->num_atoms;

        printf("attempting to exchange replicas.\n");
        for (k = gsl_rng_uniform_int(self->rng, 2);
             k <= self->num_replicas - 2;
             k += 2)
        {
                replicas_exchange(self, k);
        }

        for (size_t s = 0; s < num_iters*iters_per_cycle; s++)
                #pragma omp parallel for private(k)
                for (k = 0; k < self->num_replicas; k++)
                        simulation_next_iteration(self->replica[k]);
}

void replicas_exchange(struct replicas *self, size_t k)
{
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
                printf("swapping replicas %u and %u.\n", k, k+1);
                struct protein *x = s1->protein;
                s1->protein = s2->protein;
                s1->energy = U2;
                s2->protein = x;
                s2->energy = U1;
                ++self->exchanges[k];
        }

        ++self->total[k];

        fflush(stdout);
}



static inline bool simulation_has_converged(const struct simulation *s,
                                            double orig_energy)
{
        return gsl_fcmp(orig_energy, s->energy, 1e-3) == 0;
}

bool replicas_have_converged(const struct replicas *self)
{
        bool status = false;

        for (size_t s = 0; s < self->num_replicas; s++)
                status |= simulation_has_converged(self->replica[s],
                                                   self->orig_energy);

        return status;
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
