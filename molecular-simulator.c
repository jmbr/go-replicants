#include "molecular-simulator.h"

static bool setup_only = false;
static bool simulate_only = false;
static bool resume_simulation = false;

static struct option cmd_options[] = {
	{"resume", no_argument, (int *) &resume_simulation, true},
        {"temperature", required_argument, NULL, 't'},
        {"dmax", required_argument, NULL, 'd'},
        {"a", required_argument, NULL, 'a'},
        {"setup-only", no_argument, (int *) &setup_only, true},
        {"simulate-only", no_argument, (int *) &simulate_only, true},
        {"help", no_argument, NULL, 'h'},
        {0, 0, 0, 0}
};

static int parse_args(int argc, char *argv[], struct simulation_options *opts);
static void print_usage(void);
static void show_progress(const struct replicas *r, size_t k);


int main(int argc, char *argv[])
{
        gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_env_setup();
        gsl_rng_set(rng, gsl_rng_default_seed); /* XXX Should be a command line option. */

        const size_t max_temperatures = 256;
        double temperatures[max_temperatures];
        struct simulation_options opts = {
                rng, 0.0, 0.0, 0, (double *) &temperatures
        };

        if (parse_args(argc, argv, &opts) == -1)
                exit(EXIT_FAILURE);

        struct protein *p = protein_read_xyz_file(argv[optind++]);
        if (p == NULL) {
                perror("protein_read_xyz_file");
                exit(EXIT_FAILURE);
        }

        if (resume_simulation && argc - optind != opts.num_replicas) {
                fprintf(stderr, "The number of temperatures does not match "
                        "the number of configuration files.\n");
                exit(EXIT_FAILURE);
        }

        struct replicas *r = new_replicas(p, &opts);

        printf("Running thermalization phase.\n");

        replicas_thermalize(r, 35000);

        if (setup_only)
                exit(EXIT_SUCCESS);
        
        printf("Running production phase.\n");
        size_t k = 1;
        /* while (replicas_have_not_converged(r)) { */
        while (true) {
                replicas_next_iteration(r);
                show_progress(r, k);
                ++k;
        }

        exit(EXIT_SUCCESS);
}



int parse_args(int argc, char *argv[], struct simulation_options *opts)
{
        while (true) {
                int c = getopt_long_only(argc, argv, "", cmd_options, 0);
                if (c == -1)
                        break;
                
                switch (c) {
                case 't':
                        opts->temperatures[opts->num_replicas++] = atof(optarg);
                        break;
                case 'd':
                        opts->d_max = atof(optarg);
                        break;
                case 'a':
                        opts->a = atof(optarg);
                        break;
                case 'h':
                        print_usage();
                        exit(EXIT_SUCCESS);
                }
        }

        if (opts->d_max <= 0.0 || opts->a <= 0.0 || opts->num_replicas == 0) {
                print_usage();
                return -1;
        }

        return 0;
}

void print_usage(void)
{
        fprintf(stderr,
                "Usage: molecular-simulator [-r] -d VALUE -a VALUE "
                "-t VALUE [-t VALUE ...] PROTEIN-FILE [CONFORMATION-FILE ...]\n");
}

void show_progress(const struct replicas *r, size_t k)
{
        if (k != 1 && k % 100 != 0)
                return;

        printf("total number of exchanges: %u\n", replicas_total_exchanges(r));
        replicas_print_info(r, stdout);

        for (size_t j = 0; j < r->num_replicas; j++) {
                struct simulation *s = r->replica[j];
                simulation_print_info(s, stdout);
                printf("progress for replica %u: %g/%g\n",
                       j, s->energy, r->orig_energy);
        }

        fflush(stdout);
}
