#include <string.h>

#include "molecular-simulator.h"


int main(int argc, char *argv[])
{
        if (argc < 3) {
                fprintf(stderr, "Usage: %s NUM-ATOMS CONFORMATIONS-FILE\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        FILE *g = popen(GNUPLOT_EXECUTABLE " -persist", "w");
        if (g == NULL) {
                fprintf(stderr, "Unable to run Gnuplot.\n");
                exit(EXIT_FAILURE);
        }

        FILE *fp = fopen(argv[argc-1], "r");
        if (fp == NULL) {
                fprintf(stderr, "Unable to open `%s'.\n", argv[argc-1]);
                exit(EXIT_FAILURE);
        }

        size_t num_atoms = (size_t) atoi(argv[1]);
        double tab[num_atoms][3];

        size_t k = 0, frame = 1;
        while (true) {
                if (k == num_atoms) {
                        k = 0;
                        struct protein *p;
                        p = new_protein(num_atoms, (double *) tab);
                        protein_plot(p, g, "%s (frame %u)", argv[argc-1], frame);
                        delete_protein(p);
                }
                
                int s = fscanf(fp, "%lg %lg %lg",
                               &tab[k][0], &tab[k][1], &tab[k][2]);
                if (s != 3) {
                        fprintf(stderr, "Invalid conformation file.\n");
                        exit(EXIT_FAILURE);
                }
                
                ++k;
                ++frame;
        }

        pclose(g);

        exit(EXIT_SUCCESS);
}
