#ifndef UTILS_H
#define UTILS_H

/*
 * Some macros useful for allocating vectors on the stack.
 * Caution: these are only meant to be used with double floats.
 */
#define block_init(size_, data_) { .size = size_, .data = data_ }
#define vector_init(block, offset, n, stride_)      \
        { .data = block.data + offset,                         \
          .size = n, .stride = stride_, .owner = 0 }

/* XXX Deprecated. */
#define declare_stack_allocated_vector(name, n)                         \
        double stack_allocated_vector_array_##name[n];                  \
        gsl_vector_view stack_allocated_vector_view_##name =		\
                gsl_vector_view_array((double *) stack_allocated_vector_array_##name, n); \
        gsl_vector *name = &stack_allocated_vector_view_##name.vector

#define declare_stack_allocated_matrix(name, m, n)                      \
        double stack_allocated_matrix_array_##name[n];                  \
        gsl_matrix_view stack_allocated_matrix_view_##name =		\
                gsl_matrix_view_array((double *) stack_allocated_matrix_array_##name, m, n); \
        gsl_matrix *name = &stack_allocated_matrix_view_##name.matrix

#if DEBUG_LEVEL > 0
# define dprintf(...)                           \
do {                                            \
        fprintf(stderr, "%s: ", __func__);      \
        fprintf(stderr, __VA_ARGS__);           \
} while(0);
#else
# define dprintf(...)
#endif // !NDEBUG

#endif // !UTILS_H
