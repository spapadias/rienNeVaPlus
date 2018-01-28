//Title: winCPU.cpp
//Description: Functions for measuring the CPU time in seconds and in machine cycles
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date: 4/4/2015
//Revision: 10/5/2016

//Windows
#ifdef _WIN32

#include <Windows.h>
//#include "./Structures.h"

double get_wall_time()
{
    LARGE_INTEGER time,freq;

    if (!QueryPerformanceFrequency(&freq))
    {
        //  Handle error
        return 0;
    }
    if (!QueryPerformanceCounter(&time))
    {
        //  Handle error
        return 0;
    }
    return (double)time.QuadPart / freq.QuadPart;
}

double get_cpu_time()
{
    FILETIME a,b,c,d;
    if (GetProcessTimes(GetCurrentProcess(),&a,&b,&c,&d) != 0)
    {
        //  Returns total user time.
        //  Can be tweaked to include kernel times as well.
        return
            (double)(d.dwLowDateTime |
            ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
    }
    else
    {
        //  Handle error
        return 0;
    }
}

//Posix/Linux
#else

#include <sys/time.h>

double get_wall_time()
{
    struct timeval time;
    if (gettimeofday(&time,NULL))
    {
        //  Handle error
        return 0;
    }
    //return (double)time.tv_sec + (double)time.tv_usec * .000001;
    return (double)(time.tv_sec * 1000.0) + (double)(time.tv_usec / 1000.0);
}

double get_cpu_time()
{
    return (double)clock() / CLOCKS_PER_SEC;
}

#endif