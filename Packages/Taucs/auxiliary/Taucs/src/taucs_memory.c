/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>

#include "taucs.h"

#ifndef OSTYPE_win32
#include <unistd.h>
#endif


/********************************************************************/
/* taucs_maximize_stacksize                                         */
/*   tries to extend the stack as much as possible, to prevent      */
/*   stack overflows in recursive functions                         */
/********************************************************************/

#ifndef OSTYPE_win32
#include <unistd.h>

/* without _XPG4_2, sys/time.h does not define timeval */
#ifdef OSTYPE_solaris
#define _XPG4_2
#endif

#include <sys/time.h>
#include <sys/resource.h>

int taucs_maximize_stacksize()
{
  struct rlimit l;
  char rlim_cur[64];
  char rlim_max[64];

  if (getrlimit(RLIMIT_STACK, &l) != 0) {
    taucs_printf("taucs_maximize_stacksize: getrlimit() failed\n");
    return -1;
  }

  if (l.rlim_cur == RLIM_INFINITY) sprintf(rlim_cur,"unlimited");
  else                             sprintf(rlim_cur,"%dk",(int) l.rlim_cur / 1024);
  if (l.rlim_max == RLIM_INFINITY) sprintf(rlim_max,"unlimited");
  else                             sprintf(rlim_max,"%dk",(int) l.rlim_max / 1024);
  taucs_printf("taucs_maximize_stacksize: current stack size %s, max is %s\n",
	       rlim_cur, rlim_max);

  if (l.rlim_cur != l.rlim_max) {

    l.rlim_cur = l.rlim_max;
    
    if (setrlimit(RLIMIT_STACK, &l) != 0) {
      taucs_printf("taucs_maximize_stacksize: setrlimit() failed\n");
      return -1;
    }
    
    if (getrlimit(RLIMIT_STACK, &l) != 0) {
      taucs_printf("taucs_maximize_stacksize: getrlimit() failed\n");
      return -1;
    }
    
    if (l.rlim_cur == RLIM_INFINITY) sprintf(rlim_cur,"unlimited");
    else                             sprintf(rlim_cur,"%dk",(int) l.rlim_cur / 1024);
    if (l.rlim_max == RLIM_INFINITY) sprintf(rlim_max,"unlimited");
    else                             sprintf(rlim_max,"%dk",(int) l.rlim_max / 1024);
    taucs_printf("taucs_maximize_stacksize: current stack size %s, max is %s\n",
		 rlim_cur, rlim_max);
  }

  return 0;
}

#else /* win32 */
int taucs_maximize_stacksize()
{
  taucs_printf("taucs_maximize_stacksize: not supported on Win32,\n");
  taucs_printf("taucs_maximize_stacksize: compile with /F[stacksize] or run\n");
  taucs_printf("taucs_maximize_stacksize: EDITBIN /STACK:[stacksize] *.EXE\n");
  return -1;
}
  
#endif


/********************************************************************/
/* taucs_system_memory_size                                         */
/*   returns size of memory reported by the operating system        */
/*   should not normally be called by the user (call _avail_)       */
/********************************************************************/

#ifdef OSTYPE_linux
#define TAUCS_SYSTEM_MEMORY_SIZE_DEFINED

double taucs_system_memory_size()
{
  FILE* f;
  double m;
  double m_sc;

  m_sc  = (double) sysconf(_SC_PAGESIZE);
  m_sc *= (double) sysconf(_SC_PHYS_PAGES);

  /* total memory is the first number in /proc/meminfo */

  f = fopen("/proc/meminfo","r");
  if (f==NULL) return m_sc;
  if (fscanf(f,"%*[a-zA-Z :\n\r]%lf",&m) != 1) return m_sc;

  if (m != m_sc) {
    taucs_printf("Warning: /proc/meminfo reports %lfMB of memory while\n",
	       m/1048576.0);
    taucs_printf("         sysconf       reports %lfMB of memory\n",
	       m_sc/1048576.0);
  }

  return m;
}
#endif

#ifdef OSTYPE_darwin
#define TAUCS_SYSTEM_MEMORY_SIZE_DEFINED

/* This is a BSD4.4 interface, so it should work on other BSD systems */

#include <sys/types.h>
#include <sys/sysctl.h>

double taucs_system_memory_size()
{
  int mib[2] = { CTL_HW, HW_PHYSMEM };
  int int_retval;
  size_t len = sizeof(int);
  
  taucs_printf("taucs_system_memory_size: calling sysctl\n");
  mib[1] = HW_PAGESIZE;
  if ( sysctl(mib,2,
	      &int_retval,&len,
	      NULL, 0)) {
    taucs_printf("taucs_system_memory_size: ERROR, sysctl failed (on darwin)\n");
    return -1.0;
  }
  taucs_printf("  sysctl pagesize %d bytes\n",int_retval);

  mib[1] = HW_PHYSMEM;
  if ( sysctl(mib,2,
	      &int_retval,&len,
	      NULL, 0)) {
    taucs_printf("taucs_system_memory_size: ERROR, sysctl failed (on darwin)\n");
    return -1.0;
  }
  taucs_printf("  sysctl physmem %d bytes\n",int_retval);

  return (double) int_retval;
}
#endif

#ifdef OSTYPE_aix
#define TAUCS_SYSTEM_MEMORY_SIZE_DEFINED
double taucs_system_memory_size()
{
  FILE* f;
  double m;

  int child_stdout[2];

  pipe(child_stdout);

  if (fork() == 0) {
    char* argv[10];
    int   i = -1;
    argv[i++] = "lsattr";
    argv[i++] = "-E";
    argv[i++] = "-F";
    argv[i++] = "value";
    argv[i++] = "-l";
    argv[i++] = "sys0";
    argv[i++] = "-a";
    argv[i++] = "realmem";
    argv[i++] = 0;
    close(child_stdout[0]);
    dup2(child_stdout[1],1); 
    execv("/usr/sbin/lsattr",argv);
    perror("System error (execv)");
    exit(0);
  } else {
    char buffer[256];
    int  nbytes;
    /* parent continues */
    close(child_stdout[1]);
    nbytes = read(child_stdout[0],buffer,256);
    close(child_stdout[0]);
    if (sscanf(buffer,"%lf",&m) != 1)
      return -1.0;
    return 1024.0 * m;
  }
}
#endif

#ifdef OSTYPE_solaris
#define TAUCS_SYSTEM_MEMORY_SIZE_DEFINED
/*#include <sys/unistd.h>*/

double taucs_system_memory_size()
{
  double m;

  m  = (double) sysconf(_SC_PAGESIZE);
  m *= (double) sysconf(_SC_PHYS_PAGES);

  return m;
}
#endif

#ifdef OSTYPE_win32
#define TAUCS_SYSTEM_MEMORY_SIZE_DEFINED
#include <windows.h>

double taucs_system_memory_size()
{
  MEMORYSTATUS ms;

  GlobalMemoryStatus(&ms);

  taucs_printf("taucs_system_memory_size: returning information from GlobalMemoryStatus\n");
  taucs_printf("  Warning: may be incorrect when the machine has more than 4 GB,\n");
  taucs_printf("  Warning: or if there are more than 2 GB but the executable was\n");
  taucs_printf("  Warning: compiled without the /LARGEADDRESSWARE liner flag.\n");

  taucs_printf("  Memory load                    = %03d%%\n",ms.dwMemoryLoad);
  taucs_printf("  Total Physical Memory          = %.0f MB\n",(double) ms.dwTotalPhys /1048576.0 );
  taucs_printf("  Available Physical Memory      = %.0f MB\n",(double) ms.dwAvailPhys /1048576.0 );
  taucs_printf("  Total Page File                = %.0f MB\n",(double) ms.dwTotalPageFile /1048576.0 );
  taucs_printf("  Available Memory in Page File  = %.0f MB\n",(double) ms.dwAvailPageFile /1048576.0 );
  taucs_printf("  Address-Space Size             = %.0f MB\n",(double) ms.dwTotalVirtual /1048576.0 );
  taucs_printf("  Address-Space Available Memory = %.0f MB\n",(double) ms.dwAvailVirtual /1048576.0 );

  return (double) ms.dwTotalPhys;
}
#endif

#ifdef OSTYPE_irix
#define TAUCS_SYSTEM_MEMORY_SIZE_DEFINED
#include <sys/unistd.h>
#include <sys/types.h>
#include <sys/sysmp.h>
#include <sys/sysinfo.h>

double taucs_system_memory_size()
{
  double m;
  struct rminfo p;

  m  = (double) sysconf(_SC_PAGESIZE);
  /*m  = (double) getpagesize();*/

  taucs_printf("***************** %.0lf ************\n",m);

  if (sysmp(MP_SAGET, MPSA_RMINFO, &p) == -1) {
    perror("sysmp");
    exit(1);
  }
  m = (double) (p.physmem);
  taucs_printf("**$$$* %.0lf\n",m);  

  return m;
}
#endif

#ifndef TAUCS_SYSTEM_MEMORY_SIZE_DEFINED
#define TAUCS_SYSTEM_MEMORY_SIZE_DEFINED
double taucs_system_memory_size()
{ 
  taucs_printf("Warning: cannot automatically determine main memory size\n");
  taucs_printf("         for this platform\n");
  return -1.0; 
}
#endif


/********************************************************************/
/* taucs_avail_memory_size                                            */
/*   returns size of memory available for allocation                */
/********************************************************************/

double taucs_available_memory_size() 
{
  double m_sys;
  double m,m_low,m_high,m_tol;
  char*  p;

  m_sys = taucs_system_memory_size();
  
  /* malloc test */

  m = 1048576.0;

  while ( (p=(char*) taucs_malloc( (size_t) (m*2.0) )) != NULL ) {
    taucs_free(p);
    m = m*2.0;
  }

  m_low  = m;
  m_high = 2.0*m;
  m_tol  = m / 128.0;

  while ( m_high - m_low > m_tol ) {
    m = m_low + ( (m_high-m_low)/2.0 );
    taucs_printf("taucs_avail_memory_size: [%.0lf %.0lf %.0lf]\n",
	       m_low  / 1048576.0,
	       m      / 1048576.0,
	       m_high / 1048576.0);
    if ( (p=(char*) taucs_malloc( (size_t) m )) != NULL ) 
      m_low = m;
    else 
      m_high = m;
    taucs_free(p);
  }

  m = m_low;

  taucs_printf("taucs_avail_memory_size: malloc test=%.0lf MB sys test=%.0lf MB\n",
	     m / 1048576.0,
	     m_sys / 1048576.0
	     );

  /* if m_sys is meaningful, then we limit m by 0.75*m_sys */

  if (m_sys > 0) {
    m_sys = floor(0.75 * m_sys); 
    if (m_sys < m) m = m_sys;
  }

  return m;
}

