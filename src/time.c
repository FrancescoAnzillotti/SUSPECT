#include "suspect.h"

/* time measurement. */

double now_seconds()  //millisecons
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    // Convert seconds to ms, and nanoseconds to ms, then add them
    return (ts.tv_sec * 1000.0) + (ts.tv_nsec / 1e6);
}
