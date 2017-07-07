#pragma once
// Performance profiling macros and utilities

#define ENABLE_PERF_PROF 1

#if ENABLE_PERF_PROF

#include <chrono>

#define DECLARE_PERF_PROF_COUNTER(counter) long long counter = 0LL
#define RESET_PERF_PROF_COUNTER(counter) counter = 0LL
#define PERF_PROF_COUNTER_NS(counter) (counter)
#define START_PERF_PROF(name) auto name##_start_time = std::chrono::high_resolution_clock::now()
#define END_PERF_PROF(name) auto name##_end_time = std::chrono::high_resolution_clock::now()
#define GET_PERF_PROF_DURATION(name) get_perf_prof_duration(name##_start_time, name##_end_time)
#define ACCUM_PERF_PROF_DURATION(name, counter) counter += get_perf_prof_duration(name##_start_time, name##_end_time)

static inline long long get_perf_prof_duration(std::chrono::time_point<std::chrono::high_resolution_clock>& start_time, std::chrono::time_point<std::chrono::high_resolution_clock>& end_time)
{
	return std::chrono::duration_cast<std::chrono::nanoseconds>(end_time-start_time).count();
}

#else // if ENABLE_PERF_PROF

#define DECLARE_PERF_PROF_COUNTER(counter)
#define RESET_PERF_PROF_COUNTER(counter)
#define PERF_PROF_COUNTER_NS (-1LL)
#define START_PERF_PROF(name)
#define END_PERF_PROF(name)
#define GET_PERF_PROF_DURATION(name) (-1LL)
#define ACCUM_PERF_PROF_DURATION(name, counter)

#endif // if ENABLE_PERF_PROF
