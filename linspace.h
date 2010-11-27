#ifndef LINSPACE_H
#define LINSPACE_H

static inline void linspace(double base, double limit, size_t N, double *values)
{
        assert(base <= limit);
        assert(N > 0);
        assert(values != NULL);

        const double step = (limit - base)/(double) (N-1);

        for (size_t k = 0; k < N; k++)
                values[k] = base + k*step;
}

#endif // !LINSPACE_H
