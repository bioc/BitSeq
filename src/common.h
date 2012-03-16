#ifndef COMMON_H
#define COMMON_H

#ifdef BIOC_BUILD
#include <R.h>
#endif

//#define Sof(x) (long)x.size()
//#define FOR(x,y,n) for(x=y;x<n;x++)
//#define FR(x,n) FOR(x,0,n)
//#define FF first
//#define SS second

#ifndef BIOC_BUILD
#define message(...) printf(__VA_ARGS__)
#define errmessage(...) fprintf(stderr, __VA_ARGS__)
#define warning(x, ...) fprintf(stderr, strcat("WARNING:", x), __VA_ARGS__)
#define error(x, ...) fprintf(stderr, strcat("ERROR:", x), __VA_ARGS__)
#else
#define message(...) Rprintf(__VA_ARGS__)
#define errmessage(...) REprintf(__VA_ARGS__)
#endif

void buildTime(char *argv0);

bool progressLog(long cur,long outOf, long steps = 10);

#endif
