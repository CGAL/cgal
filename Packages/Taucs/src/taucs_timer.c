/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/

#include "taucs.h"

#ifndef TAUCS_CONFIG_TIMING
double taucs_wtime() { return 0.0; }
double taucs_ctime() { return 0.0; }
#else

#ifdef OSTYPE_win32
#define TAUCS_TIMER

#include <windows.h>

double taucs_wtime() {
  double wtime;

  SYSTEMTIME systime;

  GetSystemTime(&systime);

  wtime = 0.0
    +  0.001 * (double) systime.wMilliseconds
    +    1.0 * (double) systime.wSecond
    +   60.0 * (double) systime.wMinute
    + 3600.0 * (double) systime.wHour;

  return wtime; 
}
double taucs_ctime() { 
  double ctime;
  FILETIME creationt,exitt,kernel,user;
  ULARGE_INTEGER ukernel,uuser;
  HANDLE self;

  self = GetCurrentProcess();

  if (!GetProcessTimes(self,
		       &creationt,
		       &exitt,
		       &kernel,
		       &user))
    return 0.0;

  ukernel.LowPart  = kernel.dwLowDateTime;
  ukernel.HighPart = kernel.dwHighDateTime;
  uuser.LowPart  = user.dwLowDateTime;
  uuser.HighPart = user.dwHighDateTime;
  
  ctime = ((double) (signed __int64) ukernel.QuadPart / 10000000.0)
        + ((double) (signed __int64) uuser.QuadPart   / 10000000.0);

  CloseHandle( self );

  return ctime;
}
#endif /* win32 */

/*#include "taucs.h"*/

/* Return time in nanoseconds */
 
#if 0
#ifdef OSTYPE_linux_not_reliable
#define TAUCS_TIMER

#include <stdio.h>                                                 
#include <unistd.h>
#include <sys/types.h>                                                 
#include <sys/timeb.h>                                                 


/* p5tsc.h -- functions to use Pentium cycle counter for timing of events.
   Christian Kurmann <kurmann@inf.ethz.ch>
   based on Brad Karp, David Mazieres's p5cnt package from Harvard. */

typedef unsigned long uint32;

/* Cycle Counter */

/*
 *  Write <hi>:<lo> into MSR number <msr>.
 */

__inline__ void
wrmsr (const uint32 msr, uint32 hi, uint32 lo)
{
  __asm __volatile (
		    /*
		      "movl %0, %%ecx         # MSR to be written
		      movl %1, %%edx          # High order 32 bits
		      movl %2, %%eax          # Low order 32 bits
		      .byte 0xf; .byte 0x30   # WRMSR instruction"
		    */
        "movl %0, %%ecx movl %1, %%edx movl %2, %%eax .byte 0xf; .byte 0x30 # WRMSR instr"
    : : "g" (msr), "g" (hi), "g" (lo) : "eax", "ecx", "edx");
}

/* macro for clearing tsc */
#define cltsc wrmsr((uint32) 0x10, (uint32) (0), (uint32) (0))

/*
 *  Read 64 bit time stamp counter.  Put the high 32 bits in
 *  <*hi>, and the lower 32 bits in <*lo>.
 */
__inline__ void
rdtsc (uint32 *hi, uint32 *lo)
{
  __asm __volatile (
		    /*
		      ".byte 0xf; .byte 0x31  # RDTSC instruction
		      movl %%edx, %0          # High order 32 bits
		      movl %%eax, %1          # Low order 32 bits"
		    */
      ".byte 0xf; .byte 0x31 movl %%edx, %0 movl %%eax, %1 # RDTSC instruction"
    : "=g" (*hi), "=g" (*lo) :: "eax", "edx");
}

/* Performance Monitor Counters */

__inline__ void
spmc (uint32 ecx)
{
  __asm __volatile (
        "movl %0, %%ecx         # select counter "
    : : "g" (ecx) : "ecx");
}

/*
 *  Read 64 bit Performance Monitor Counter.  Put the high 32 bits in
 *  <*hi>, and the lower 32 bits in <*lo>.
 */
__inline__ void
rdpmc (uint32 *hi, uint32 *lo)
{
  __asm __volatile (
		    /*
		      ".byte 0xf; .byte 0x33  # RDPMC instruction
		      movl %%edx, %0          # High order 32 bits
		      movl %%eax, %1          # Low order 32 bits"
		    */
    ".byte 0xf; .byte 0x33 movl %%edx, %0 movl %%eax, %1 # RDPMC instruction"
    : "=g" (*hi), "=g" (*lo) :: "eax", "edx");
}

/* 64 bit subtract t1-t0 (result 32 bit integer) */
int subtract64(uint32 hi0, uint32 lo0, uint32 hi1, uint32 lo1 )
{
  uint32 hir, lor;

  hir = (hi1 - hi0);
  lor = (lo1 - lo0);
  if (lo1 < lo0) hir -= 1;
  return (hir > 0 ? 0:lor);
}

double timer()
{
  uint32        hi, lo;
  static uint32 hi0, lo0;
  static uint32 hi1, lo1;
    uint32 hir, lor;
  static double loticks_per_ns;
  static double hiticks_per_ns;
  static int    first_time    = 1;

  if (first_time) {
    struct timeb T;
    static time_t start_time, time_diff;
    static time_t start_mill, mill_diff;
    int    rc;
    double dt;
    /*uint32 hi0, lo0, hi1, lo1;*/
    double ticks;

    first_time = 0;

    rc = ftime( &T );
    rdtsc(&hi0,&lo0);
    start_time = T.time;
    start_mill = T.millitm;

    sleep(1);

    rc = ftime( &T );
    rdtsc(&hi1,&lo1);
    time_diff = T.time - start_time;
    mill_diff = T.millitm - start_mill; 

    dt = (1e9) * ((double) time_diff) + (1e6) * ((double) mill_diff);

    hir = (hi1 - hi0);
    lor = (lo1 - lo0);
    if (lo1 < lo0) hir -= 1;
    ticks = (double) lor + (double) hir * 4294967296.0;

    loticks_per_ns = ticks/dt;
    hiticks_per_ns = ticks/(dt / 4294967296.0);

    hir = (hi1 - hi0);
    loticks_per_ns = (4294967296.0 * (double) hir + (double) lo1 - (double) lo0) / dt;

    log_printf ("timer: hi0 %u lo0 %u hi1 %u lo1 %u\n",hi0,lo0,hi1,lo1);
    log_printf ("timer: lo/ns %lg hi/ns %lg\n",loticks_per_ns,hiticks_per_ns);

    rdtsc(&hi0,&lo0);
  }

  rdtsc(&hi1,&lo1);

  hir = (hi1 - hi0);
  return ( loticks_per_ns * ( 4294967296.0 * (double) hir + (double) lo1 ) );

  /*
  rdtsc(&hi,&lo);
  return ( ((double) hi * hiticks_per_ns) + ((double) lo * loticks_per_ns) );
  */
}

#endif /* OSTYPE_linux_not_reliable */
#endif /* if 0 */

/*********************************************************/
/*                                                       */
/*********************************************************/

#ifndef OSTYPE_win32
#define TAUCS_TIMER

#include <stdio.h>                                                 
#include <unistd.h>

/*
  #ifndef OSTYPE_solaris
  #include <sys/time.h>
  #include <sys/resource.h>
  #endif
*/
#ifdef OSTYPE_solaris
#define _XPG4_2
#endif
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>                                                 
#include <sys/timeb.h>                                                 

double taucs_wtime()
{
  struct timeb T;
  /*static int first_time    = 1;*/
  /*  static time_t start_time, time_diff;
      static time_t start_mill, mill_diff;
  */

  static time_t time_diff;
  static time_t mill_diff;
  /*int    rc;*/
  double dt;
  
  (void) ftime( &T );
  /*
  if (first_time) {
    first_time = 0;
    start_time = T.time;
    start_mill = T.millitm;
  }

  time_diff = T.time - start_time;
  mill_diff = T.millitm - start_mill; 
  */
  time_diff = T.time;
  mill_diff = T.millitm;

  dt = ((double) time_diff) + (1e-3) * ((double) mill_diff);

  return dt;
}

double taucs_ctime()
{
  /*
  #ifdef OSTYPE_solaris
  return 0.0;
  #else
  */
  struct rusage a;
  
  getrusage(RUSAGE_SELF,&a);

  return (double) 
    (double) ((a.ru_utime).tv_sec +(a.ru_stime).tv_sec ) +
    (double) ((a.ru_utime).tv_usec+(a.ru_stime).tv_usec) * 1.0e-6;
  /*#endif*/
}

#endif /* not win32 */

/*********************************************************/
/*                                                       */
/*********************************************************/

#ifndef TAUCS_TIMER
#define TAUCS_TIMER
double taucs_wtime() { return 0.0; }
double taucs_ctime() { return 0.0; }

#endif

#if 0
void cpu_time_from_last(char s[])
{
  struct rusage a;
  struct timeb T;
  static int first_time    = 1;
  static time_t start_time, time_diff;
  static time_t start_mill, mill_diff;
  /*int    rc;*/
  double dt,cpu_t;
  static double last_cpu_t;

  (void) ftime( &T );

  if (first_time) {
    first_time = 0;
    start_time = T.time;
    start_mill = T.millitm;
    getrusage(RUSAGE_SELF,&a);
    
    last_cpu_t = (a.ru_utime).tv_sec+(a.ru_stime).tv_sec+
      ((a.ru_utime).tv_usec+(a.ru_stime).tv_usec)/1000000.0; 
      

    taucs_printf("starting timer\n");
  }
  else
    {
      time_diff = T.time - start_time;
      mill_diff = T.millitm - start_mill; 
      
      dt = ((double) time_diff) + ((double) mill_diff)/1000.;

      start_time = T.time;
      start_mill = T.millitm;
      
      getrusage(RUSAGE_SELF,&a);

      cpu_t = (a.ru_utime).tv_sec+(a.ru_stime).tv_sec+
	((a.ru_utime).tv_usec+(a.ru_stime).tv_usec)/1000000.0; 
      
      taucs_printf("%s - CPU Time from last time : %lf\n", s, cpu_t-last_cpu_t);
      taucs_printf("%s - Wall Clock Time from last time : %lf\n",s,dt);
      
      last_cpu_t = cpu_t;
    }
}

void cpu_time_from_start(char s[])
{
  struct rusage a;
  struct timeb T;
  static int first_time    = 1;
  static time_t start_time, time_diff;
  static time_t start_mill, mill_diff;
  /*int    rc;*/
  double dt,cpu_t;
  static double start_cpu_t;
  
  (void) ftime( &T );

  if (first_time) {
    first_time = 0;
    start_time = T.time;
    start_mill = T.millitm;
    getrusage(RUSAGE_SELF,&a);

    start_cpu_t = (a.ru_utime).tv_sec+(a.ru_stime).tv_sec+
      ((a.ru_utime).tv_usec+(a.ru_stime).tv_usec)/1000000.0; 
      

    taucs_printf("starting timer\n");
  }
  else
    {
      time_diff = T.time - start_time;
      mill_diff = T.millitm - start_mill; 
      
      dt = ((double) time_diff) + ((double) mill_diff)/1000.;
      
      getrusage(RUSAGE_SELF,&a);

      cpu_t = (a.ru_utime).tv_sec+(a.ru_stime).tv_sec+
	((a.ru_utime).tv_usec+(a.ru_stime).tv_usec)/1000000.0; 
      
      taucs_printf("%s - CPU Time from beginning : %lf\n",s,cpu_t-start_cpu_t);
      taucs_printf("%s - Wall Clock Time from beginning : %lf\n",s,dt);
      
    }
}
#endif /* if 0 */

#endif /* TAUCS_CONFIG_TIMING */

/*********************************************************/
/*                                                       */
/*********************************************************/
