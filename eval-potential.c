#include "molecular-simulator.h"


static void print_usage(void);


int main(int argc, char *argv[])
{
        set_prog_name("eval-potential");

        char *reference = NULL;
        double d_max = 0.0, a = 0.0;

        while (true) {
                struct option cmd_options[] = {
                        {"reference", required_argument, NULL, 'r'},
                        {"dmax", required_argument, NULL, 'd'},
                        {"a", required_argument, NULL, 'a'},
                        {"help", no_argument, NULL, 'h'},
                        {0, 0, 0, 0}
                };
                
                int c = getopt_long_only(argc, argv, "", cmd_options, 0);
                if (c == -1)
                        break;

                switch (c) {
                case 'r':
                        reference = optarg;
                        break;
                case 'd':
                        d_max = atof(optarg);
                        break;
                case 'a':
                        a = atof(optarg);
                        break;
                case 'h':
                        print_usage();
                        exit(EXIT_SUCCESS);
                }
        }

        if (reference == NULL || d_max <= 0.0 || a <= 0.0 || optind == argc) {
                print_usage();
                exit(EXIT_FAILURE);
        }
        
        struct protein *p1, *p2;
        struct contact_map *m1, *m2;

        p1 = protein_read_xyz_file(reference);
        if (p1 == NULL)
                die_printf("Unable to open `%s'.\n", reference);


        FILE *f = fopen(argv[optind], "r");
        if (f == NULL)
                die_printf("Unable to open `%s'.\n", argv[optind]);
        p2 = protein_read_latest_xyz(f);
        fclose(f);

        m1 = new_contact_map(p1, d_max);
        m2 = new_contact_map(p2, d_max);
        if (m1 == NULL || m2 == NULL)
                die("Unable to compute contact maps.");

        contact_map_diff(m1, m2);

        printf("%f\n", potential(p2, m1, a));

        exit(EXIT_SUCCESS);
}

void print_usage(void)
{
        printf("Usage: %s --reference FILE --dmax VALUE "
               "--a VALUE FILE\n", get_prog_name());
}
