#undef NDEBUG
#include "molecular-simulator.h"


static bool plot_results = false;


static void test_initialization_and_finalization(void);
static void test_potential_energy(void);


int main(int argc, char *argv[])
{
        int opt;
        
        while ((opt = getopt(argc, argv, "gh")) != -1) {
                switch (opt) {
                case 'g':
                        plot_results = true;
                        break;
                case 'h':
                        fprintf(stderr, "Usage: %s [-g]\n", argv[0]);
                        exit(EXIT_SUCCESS);
                }
        }

        test_initialization_and_finalization();
        test_potential_energy();

        exit(EXIT_SUCCESS);
}


void test_initialization_and_finalization(void)
{
        struct protein *p = new_protein_2gb1();
        assert(p != NULL);

        struct contact_map *c1 = new_contact_map(p, 10);
        struct contact_map *c2 = new_contact_map(p, 7.5);
        assert(c1 != NULL && c2 != NULL);
        
        if (plot_results) {
                printf("The following plots should be compared against those appearing"
                       "in Figure 3.2 of Lidia Prieto's Ph.D. thesis (p. 60).\n");

                FILE *g = popen("gnuplot --persist", "w");
                assert(g != NULL);

                contact_map_plot(c1, g, "");
                getchar();
                contact_map_plot(c2, g, "");

                pclose(g);
        }

        delete_contact_map(c1);
        delete_contact_map(c2);
        delete_protein(p);
}

void test_potential_energy(void)
{
        struct protein *p1 = new_protein_1pgb();
        struct protein *p2 = new_protein_2gb1();
        assert(p1 != NULL && p2 != NULL);
        
        const double d_max = 10;
        struct contact_map *c1 = new_contact_map(p1, d_max);
        struct contact_map *c2 = new_contact_map(p2, d_max);
        assert(c1 != NULL && c2 != NULL);

        double p_1pgb = potential(p1, c1, 0.9);
        double p_2gb1 = potential(p2, c2, 0.9);

        printf("potential(1pgb) = %g\n", p_1pgb);
        printf("potential(2gb1) = %g\n", p_2gb1);

        assert(gsl_fcmp(p_1pgb, -460.0, 1e-15) == 0);
        assert(gsl_fcmp(p_2gb1, -460.0, 1e-15) == 0);

        delete_contact_map(c1);
        delete_contact_map(c2);
        delete_protein(p1);
        delete_protein(p2);
}
