#ifndef __POLILOG_H_
#define __POLILOG_H_

#ifdef __cplusplus
extern "C"
{
#endif

#include "PoLiMEr.h"

#define RED     "\x1b[31m"
#define GREEN   "\x1b[32m"
#define BLUE    "\x1b[34m"
#define YELLOW  "\x1b[33m"
#define CYAN    "\x1b[36m"
#define MAGENTA "\x1b[35m"
#define RESET   "\x1b[0m"

typedef enum levels {
    CRITICAL = 10,
    ERROR    = 20,
    WARNING  = 30,
    INFO     = 40,
    DEBUG    = 50,
    TRACE    = 60
} poliloglevel_t;

void poli_setloglevel(poliloglevel_t level);

void poli_log(poliloglevel_t level, struct monitor_t * monitor, char *format_str, ...);

#ifdef __cplusplus
}
#endif

#endif