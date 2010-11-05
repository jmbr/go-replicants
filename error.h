#ifndef ERROR_H
#define ERROR_H
/**
 * @file error.h
 * @brief Error reporting functions.
 */

extern void info_printf(const char *format, ...);

extern void error_printf(const char *format, ...);

extern void fatal_error(const char *format, ...) __attribute__ ((noreturn));
#endif // !ERROR_H
