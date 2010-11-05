/**
 * @file error.c
 * @brief Error reporting functions.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "error.h"


void info_printf(const char *format, ...)
{
        fprintf(stdout, "Info: ");

	va_list ap;
	va_start(ap, format);
	vfprintf(stdout, format, ap);
        fflush(stdout);
	va_end(ap);
}


void error_printf(const char *format, ...)
{
        fprintf(stderr, "Error: ");

	va_list ap;
	va_start(ap, format);
	vfprintf(stderr, format, ap);
        fflush(stderr);
	va_end(ap);
}

void fatal_error(const char *format, ...)
{
        fprintf(stderr, "Fatal error: ");

	va_list ap;
	va_start(ap, format);
	vfprintf(stderr, format, ap);
        fflush(stderr);
	va_end(ap);

	exit(EXIT_FAILURE);
}
