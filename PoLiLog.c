#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "PoLiLog.h"

#ifdef _DEBUG
static int poliloglevel = DEBUG;
#else
#ifdef _TRACE
static int poliloglevel = TRACE;
#else
static int poliloglevel = WARNING;
#endif
#endif

void poli_setloglevel(poliloglevel_t level)
{
    poliloglevel = level;
}

void poli_log(poliloglevel_t level, struct monitor_t * monitor, char *format_str, ...)
{
    char *levelstr, *color;

    va_list argptr;

    if (level > poliloglevel)
        return;

    switch (level)
    {
        case CRITICAL:
            levelstr = "CRITICAL";
            color = MAGENTA;
            break;
        case ERROR:
            levelstr = "ERROR";
            color = RED;
            break;
        case WARNING:
            levelstr = "WARNING";
            color = YELLOW;
            break;
        case INFO:
            levelstr = "INFO";
            color = GREEN;
            break;
        case DEBUG:
            levelstr = "DEBUG";
            color = CYAN;
            break;
        case TRACE:
            levelstr = "TRACE";
            color = BLUE;
            break;
    }

    if (level == CRITICAL || level == ERROR)
    {
        fprintf(stderr, "%s [%s]" RESET, color, levelstr);
        if (monitor != NULL)
            fprintf(stderr, "Node: %s, Rank: %d --> ", monitor->my_host, monitor->world_rank);
        va_start(argptr, format_str);
        vfprintf(stderr, format_str, argptr);
        fprintf(stderr, "\n");
        va_end(argptr);
    }
    else
    {
        printf("%s [%s]" RESET, color, levelstr);
        if (monitor != NULL)
            printf(" Node: %s, Rank: %d --> ", monitor->my_host, monitor->world_rank);
        va_start(argptr, format_str);
        vprintf(format_str, argptr);
        printf("\n");
        va_end(argptr);
    }
}