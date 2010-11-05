/**
 * @file xalloc.c
 * @brief Error-checking wrappers for memory allocation functions.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <strings.h>
#include <errno.h>

#include "error.h"
#include "xalloc.h"


void *xmalloc(size_t size)
{
	void *ptr;

	if ((ptr = malloc(size)) == NULL)
		fatal_error("malloc: %s", strerror(errno));

	memset(ptr, 0, size);

	return ptr;
}

void *xcalloc(size_t nmemb, size_t size)
{
	void *ptr;

	if ((ptr = calloc(nmemb, size)) == NULL)
		fatal_error("calloc: %s", strerror(errno));

	return ptr;
}

void *xrealloc(void *ptr, size_t size)
{
	void *p;

	if ((p = realloc(ptr, size)) == NULL)
		fatal_error("realloc: %s", strerror(errno));

	return p;
}

void xfree(void *ptr)
{
        if (ptr == NULL)
                fatal_error("free: null pointer.");
        free(ptr);
}


/* char *xstrdup(const char *s) */
/* { */
/* 	char *p; */

/* 	if ((p = strdup(s)) == NULL) */
/* 		fatal_error("strdup: %s", strerror(errno)); */

/* 	return p; */
/* } */
