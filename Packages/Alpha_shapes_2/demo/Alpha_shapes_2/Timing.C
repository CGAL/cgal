#include <stdio.h>
#include        "Timing.h"

#ifndef CLK_TCK
#define CLK_TCK 100
#endif

extern "C" {

/*
 * Structure returned by times()
 */
struct tms {
	clock_t	tms_utime;		/* user time */
	clock_t	tms_stime;		/* system time */
	clock_t	tms_cutime;		/* user time, children */
	clock_t	tms_cstime;		/* system time, children */
};

#if !defined(_KERNEL)
#if defined(__STDC__)
  clock_t times(struct tms *);
#else
  clock_t times();
#endif
#endif

}
//---------------------------- Timing functions -----------------------------

struct tms starttime, endtime;

void start_timing(tms  &starttime)
{
  times(&starttime);
}

void start_timing()
{
  times(&starttime);
}

// CLK_TCK  clock ticks per second 
void end_timing(int flag)
{
  times(&endtime);
  if (flag != 0)
    printf("TIMING %d:  %ld ticks  %.3f secs \n", flag,
	   (endtime.tms_utime - starttime.tms_utime),
	   (float) (endtime.tms_utime - starttime.tms_utime) / CLK_TCK);
}


void end_timing(int flag, const tms &starttime)
{
  times(&endtime);
  if (flag != 0)
    printf("TIMING %d:  %ld ticks  %.3f secs \n", flag,
	   (endtime.tms_utime - starttime.tms_utime),
	   (float) (endtime.tms_utime - starttime.tms_utime) / CLK_TCK);
// does not this depend an the frequency of the machine
// f = 1/T
}

void cumulate_timing(long &cumulated_utime, const tms &starttime)
{
  times(&endtime);
  cumulated_utime +=(endtime.tms_utime - starttime.tms_utime);
}

void init_cumulate_timing(long &cumulated_utime)
{
  cumulated_utime=0;
}

void print_cumulate_timing(const long &cumulated_utime, int flag)
{
  if (flag != 0)
    printf("TIMING %d:  %ld ticks  %.3f secs \n", flag,
	   cumulated_utime,
	   (float) cumulated_utime / CLK_TCK);
}
