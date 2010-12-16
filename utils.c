#include "molecular-simulator.h"


static const char *prog_name;


void set_prog_name(const char *name)
{
        assert(name != NULL);
        prog_name = name;
}

const char *get_prog_name(void)
{
        return prog_name;
}


void die(const char *message)
{
        fprintf(stderr, "%s: %s\n", prog_name, message);
        exit(EXIT_FAILURE);
}

void die_errno(const char *func_name)
{
        fprintf(stderr, "%s: %s: %s\n", prog_name, func_name, strerror(errno));
        exit(EXIT_FAILURE);
}

void die_printf(const char *format, ...)
{
        fprintf(stderr, "%s: ", prog_name);

        va_list ap;
        va_start(ap, format);
        vfprintf(stderr, format, ap);
        va_end(ap);

        exit(EXIT_FAILURE);
}
