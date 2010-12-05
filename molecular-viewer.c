#include "molecular-simulator.h"


int main(int argc, char *argv[])
{
        if (argc < 2) {
                fprintf(stderr, "Usage: %s [-d VALUE] FILE\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        double d_max = 0.0;
        if (strcmp(argv[1], "-d") == 0)
                d_max = atof(argv[2]);
        
        struct protein *p = protein_read_xyz_file(argv[argc-1]);
        if (p == NULL) {
                fprintf(stderr, "Unable to read `%s'.\n", argv[argc-1]);
                exit(EXIT_FAILURE);
        }

        FILE *g1 = popen(GNUPLOT_EXECUTABLE " -persist", "w");
        if (g1 == NULL) {
                fprintf(stderr, "Unable to run Gnuplot.\n");
                exit(EXIT_FAILURE);
        }

        protein_plot(p, g1, "%s", argv[argc-1]);

        pclose(g1);

        if (d_max > 0.0) {
                struct contact_map *m = new_contact_map(p, d_max);
                FILE *g2 = popen(GNUPLOT_EXECUTABLE " -persist", "w");
                if (g2 == NULL) {
                        fprintf(stderr, "Unable to run Gnuplot.\n");
                        exit(EXIT_FAILURE);
                }

                contact_map_plot(m, g2,
                                 "%s (native contacts: %u, pot. energy: %f)",
                                 argv[argc-1], contact_map_get_num_contacts(m),
                                 potential(p, m, 0.1));

                pclose(g2);
        }

        exit(EXIT_SUCCESS);
}
