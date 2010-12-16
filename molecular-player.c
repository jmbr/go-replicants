#include "molecular-simulator.h"


int main(int argc, char *argv[])
{
        set_prog_name("molecular-player");

        if (argc < 2) {
                fprintf(stderr, "Usage: %s FILE\n", get_prog_name());
                exit(EXIT_FAILURE);
        }

        FILE *g = popen(GNUPLOT_EXECUTABLE " -persist", "w");
        if (g == NULL) die("Unable to run Gnuplot.");

        FILE *fp = fopen(argv[1], "r");
        if (fp == NULL) die_errno("fopen");

        size_t frame = 1;
        while (true) {
                struct protein *p = protein_read_xyz(fp);
                if (feof(fp))
                        break;
                if (p == NULL) die_errno("protein_read_xyz");
                protein_plot(p, g, false, "%s (frame %u)", argv[1], frame);
                delete_protein(p);
                ++frame;
        }

        pclose(g);
        exit(EXIT_SUCCESS);
}
