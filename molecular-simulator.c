#include "molecular-simulator.h"


static void print_usage(void);
static void show_progress(const struct replicas *r, size_t k);


int main(int argc, char *argv[])
{
        set_prog_name("molecular-simulator");

        bool setup_only = false, simulate_only = false, resume = false;
        gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_env_setup();
        gsl_rng_set(rng, gsl_rng_default_seed);
        const size_t max_temperatures = 256;
        double temperatures[max_temperatures];
        struct simulation_options opts = {
                .rng = rng, .d_max = 0.0, .a = 0.0,
                .num_replicas = 0, .temperatures = (double *) &temperatures
        };

        while (true) {
                struct option long_options[] = {
                        {"resume", no_argument, (int *) &resume, true},
                        {"temperature", required_argument, NULL, 't'},
                        {"dmax", required_argument, NULL, 'd'},
                        {"a", required_argument, NULL, 'a'},
                        {"setup-only", no_argument, (int *) &setup_only, true},

                        {"simulate-only", no_argument, (int *) &simulate_only, true},
                        {"help", no_argument, NULL, 'h'},
                        {0, 0, 0, 0}
                };

                int c = getopt_long_only(argc, argv, "", long_options, NULL);

                if (c == -1)
                        break;

                switch (c) {
                case 't':
                        temperatures[opts.num_replicas++] = atof(optarg);
                        break;
                case 'd':
                        opts.d_max = atof(optarg);
                        break;
                case 'a':
                        opts.a = atof(optarg);
                        break;
                case 'h':
                        print_usage();
                        exit(EXIT_SUCCESS);
                }
        }

        if (opts.d_max <= 0.0 || opts.a <= 0.0 || opts.num_replicas == 0
            || (setup_only && simulate_only))
        {
                print_usage();
                exit(EXIT_FAILURE);
        }

        if (resume && (argc - optind - 1 != (int) opts.num_replicas))
                die("The number of temperatures does not match "
                    "the number of configuration files.");

        struct replicas *r;
        r = new_replicas(protein_read_xyz_file(argv[optind++]), &opts);
        if (r == NULL)
                die_printf("Unable to set up replicas (%s).\n", strerror(errno));

        if (!resume) {
                replicas_first_iteration(r);
        } else {
                struct protein *X[opts.num_replicas];

                for (size_t k = 0; k < opts.num_replicas; k++) {
                        printf("Reading `%s'.\n", argv[optind + (int) k]);
                        char *name = argv[optind + (int) k];
                        FILE *f = fopen(name, "r");
                        if (f == NULL)
                                die_printf("Unable to open `%s'.\n", name);
                        X[k] = protein_read_latest_xyz(f);
                        fclose(f);
                }

                replicas_resume(r, (const struct protein **) X);

                for (size_t k = 0; k < opts.num_replicas; k++)
                        delete_protein(X[k]);
        }
        
        if (!simulate_only) {
                printf("Running thermalization phase.\n");
                replicas_thermalize(r, 2500000);
                if (setup_only) {
                        printf("Finished setup phase.\n");
                        exit(EXIT_SUCCESS);
                }
        }
        
        printf("Running production phase.\n");
        size_t k = 1;
        /* while (replicas_have_not_converged(r)) { */
        while (true) {
                replicas_next_iteration(r);
                show_progress(r, k);
                ++k;
        }

        delete_replicas(r);
        gsl_rng_free(rng);
        exit(EXIT_SUCCESS);
}



void print_usage(void)
{
        fprintf(stderr,
                "Usage: molecular-simulator [--resume] [--setup-only] [--simulate-only] "
                "-d VALUE -a VALUE -t VALUE [-t VALUE ...] PROTEIN-FILE "
                "[CONFORMATION-FILE ...]\n");
}

void show_progress(const struct replicas *r, size_t k)
{
        if (k != 1 && k % 100 != 0)
                return;

        printf("total number of exchanges: %zu\n", replicas_total_exchanges(r));
        replicas_print_info(r, stdout);

        for (size_t j = 0; j < r->num_replicas; j++) {
                struct simulation *s = r->replica[j];
                simulation_print_info(s, stdout);
                printf("replica %zu: %g\n", j, s->energy);
        }

        fflush(stdout);
}
