#include "molecular-simulator.h"


static struct {
        bool plot_results;
        char *resume_file;
        double temperature;
} options = { false, NULL, GSL_NEGINF };

static FILE *gp;


static int process_cmd_args(int argc, char *argv[],
                            struct protein **orig, struct protein **initial);
static void print_usage(void);
static void show_progress(struct simulation *s, size_t k);



static void plot_protein_and_sleep(const struct simulation *s)
{
        if (options.plot_results) {
                protein_plot(s->protein, s->gnuplot,
                             "Potential energy: %g.", s->energy);
                sleep(3);
        }
}

int main(int argc, char *argv[])
{
        struct protein *orig = NULL, *initial = NULL;

        if (process_cmd_args(argc, argv, &orig, &initial) == -1) {
                if (orig)
                        delete_protein(orig);
                if (initial)
                        delete_protein(initial);
                exit(EXIT_FAILURE);
        }

        gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_env_setup();
        gsl_rng_set(rng, gsl_rng_default_seed);

        struct simulation_options opts = { .d_max = 7.5, .a = 1.0 };
        struct simulation *s = new_simulation(orig, rng, options.temperature, &opts);
        if (s == NULL) {
                fprintf(stderr, "Unable to set up simulation.\n");
                exit(EXIT_FAILURE);
        }

        plot_protein_and_sleep(s);

        simulation_print_info(s, stdout);

        if ((gp = popen("gnuplot -p", "w")) == NULL) {
                fprintf(stderr, "Unable to launch Gnuplot.\n");
                exit(EXIT_FAILURE);
        }
        
        if (initial != NULL) {
                delete_protein(s->protein);
                s->protein = initial;
                s->energy = potential(s->protein, s->native_map, s->a);
        } else {
                simulation_first_iteration(s);
        }

        plot_protein_and_sleep(s);

        size_t k = 1;
        /* while (simulation_has_not_converged(s)) { */
        while (true) {
                simulation_next_iteration(s);
                show_progress(s, k);
                ++k;
        }
        
        fclose(gp);
        delete_protein(s->protein);
        delete_simulation(s);
        exit(EXIT_SUCCESS);
}



int process_cmd_args(int argc, char *argv[],
                     struct protein **orig, struct protein **initial)
{
        int opt;
        while ((opt = getopt(argc, argv, "t:r:gh")) != -1) {
                switch (opt) {
                case 't':
                        options.temperature = atof(optarg);
                        break;
                case 'r':
                        options.resume_file = optarg;
                        break;
                case 'g':
                        options.plot_results = true;
                        break;
                case 'h':
                        print_usage();
                        exit(EXIT_SUCCESS);
                }
        }

        if (optind >= argc) {
                print_usage();
                return -1;
        }

        if (options.temperature == GSL_NEGINF) {
                fprintf(stderr, "You must specify a value for the temperature.\n");
                print_usage();
                return -1;
        }

        if ((*orig = protein_read_xyz_file(argv[optind])) == NULL) {
                fprintf(stderr, "Unable to read %s.\n", argv[optind]);
                return -1;
        }

        if (options.resume_file != NULL) {
                if ((*initial = protein_read_xyz_file(options.resume_file)) == NULL) {
                        fprintf(stderr, "Unable to read %s.\n", options.resume_file);
                        return -1;
                }
        }

        return 0;
}

void print_usage(void)
{
        fprintf(stderr, "Usage: molecular-simulator [-g] [-r restart-file] source-file\n");
}



void show_progress(struct simulation *s, size_t k)
{
        if (k != 1 && k % 10000 != 0)
                return;

        protein_write_xyz_file(s->protein, "latest.xyz");

        if (options.plot_results) {
                protein_plot(s->protein, s->gnuplot,
                             "#%u (%g).  Acceptance: %2.03g",
                             k, s->energy, simulation_get_acceptance_ratio(s));

                if (k % 100000 == 0)
                        fprintf(gp,
                        "set terminal wxt noraise\n"
                        "set grid\n"
                        "plot 'energies.dat' title 'Potential energy' with lines\n");
        }
        
        printf("Progress %1.3g%% (%g/%g)                    \r",
               100*s->energy/s->orig_energy, s->energy, s->orig_energy);
}
