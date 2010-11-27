#include "molecular-simulator.h"


static inline double pairwise_potential(const struct protein *self,
                                        size_t i, size_t j,
                                        double a, double dnat);


double potential(const struct protein *p,
                 const struct contact_map *native_map,
                 double a)
{
        assert(p != NULL);
        assert(a > 0.0);

        double U = 0.0;
        
        /* XXX It is not necessary to check all contacts in each
         * iteration (just those that are finite). */

        for (size_t i = 0; i < p->num_atoms; i++) {
                for (size_t j = i; j < p->num_atoms; j++) {
                        double d_nat = contact_map_get_distance(native_map, i, j);

                        if (gsl_isinf(d_nat))
                                continue;
                        
                        U += pairwise_potential(p, i, j, a, d_nat);
                }
        }

        return U;
}

static inline double pairwise_potential(const struct protein *p,
                                        size_t i, size_t j,
                                        double a, double dnat)
{
        const double a2 = gsl_pow_2(a);
        const double r = protein_signum(p, i, j)*protein_distance(p, i, j);

        return (fabs(r - dnat) < a) ? 1/a2*(gsl_pow_2(r - dnat) - a2) : 0;
}
