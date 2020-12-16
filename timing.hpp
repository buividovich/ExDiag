#ifndef _TIMING_HPP_
#define _TIMING_HPP_

#include <omp.h>

#define TIMING_START  {a_time = omp_get_wtime();};
#define TIMING_FINISH {a_time = omp_get_wtime() - a_time;};
#define TIMING_END    {a_time = omp_get_wtime() - a_time;};
#define TIMING_STOP   {a_time = omp_get_wtime() - a_time;};

#endif
