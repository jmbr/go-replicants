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
        size_t i, j;
        const size_t N = p->num_atoms;
#pragma omp parallel for private(i, j) reduction(+:U)
        for (i = 0; i < N; i++) {
                for (j = i+2; j < N; j++) {
                        double d_nat = contact_map_get_distance(native_map, i, j);

                        if (d_nat == 0.0)
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
        const double r = protein_signum(p, i, j)*protein_distance(p, i, j);

        return fabs(r - dnat) < a ? -1.0 + gsl_pow_2((r - dnat)/a)
                                  :  0.0;
}
