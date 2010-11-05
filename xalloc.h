#ifndef XALLOC_H
# define XALLOC_H
/**
 * @file xalloc.c
 * @brief Error-checking wrappers for memory allocation functions.
 */

#include <sys/types.h>

extern void *xmalloc(size_t size);

extern void *xcalloc(size_t nmemb, size_t size);

extern void *xrealloc(void *ptr, size_t size);

extern void xfree(void *ptr);

extern char *xstrdup(const char *s);
#endif // !XALLOC_H
