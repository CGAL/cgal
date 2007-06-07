// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy


#ifndef CGAL_TAUCS_FIX
#define CGAL_TAUCS_FIX

// Include TAUCS main header taucs.h
#ifndef CGAL_INCLUDED_TAUCS_H
    #define CGAL_INCLUDED_TAUCS_H

    // GCC 3.x does not compile if we include complex.h within "extern C {}"
    // and complains that this header is deprecated.
    #if defined(__GNUC__)
        #undef __DEPRECATED
        #include <complex.h>
    #endif

    // taucs.h will define min/max macros if it's not already done (e.g. by Windows.h).
    #ifndef min
        #define CGAL_TAUCS_DEFINES_MIN
    #endif
    #ifndef max
        #define CGAL_TAUCS_DEFINES_MAX
    #endif

    // TAUCS is a C library
    extern "C" {
        #include <taucs.h>
    }

    // Undefine Taucs' min/max macros to avoid an error
    // with std::min()/std::max() calls in standard C++ headers.
    #ifdef CGAL_TAUCS_DEFINES_MIN
        #undef min
    #endif
    #ifdef CGAL_TAUCS_DEFINES_MAX
        #undef max
    #endif
#endif

#ifdef OSTYPE_linux
    #include <unistd.h>
    #include <math.h>
#endif


extern "C"
{


/********************************************************************/
/* FIX OF taucs_system_memory_size()                                */
/*   returns size of memory reported by the operating system        */
/*   should not normally be called by the user (call _avail_)       */
/* (see taucs_memory.c)                                             */
/********************************************************************/

#ifdef OSTYPE_linux

/* Redirect call to avoid link error */
#define taucs_system_memory_size cgal_taucs_system_memory_size

inline double taucs_system_memory_size()
{
  /* LS 2006: The original code below is buggy on Linux 2.6 */
  /*          (because /proc/meminfo format changed)        */
#if 0
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
#endif
  /* LS 2006: fix */
  double m_sc;
  m_sc  = (double) sysconf(_SC_PAGESIZE);
  m_sc *= (double) sysconf(_SC_PHYS_PAGES);
  return m_sc;
}
#endif

/* LS 2007: Cygwin implementation */
#ifdef OSTYPE_cygwin

/* Redirect call to avoid link error */
#define taucs_system_memory_size cgal_taucs_system_memory_size

inline double taucs_system_memory_size()
{
  double m_sc;
  m_sc  = (double) 4096; /* page size */
  m_sc *= (double) sysconf(_SC_PHYS_PAGES);
  return m_sc;
}
#endif



/********************************************************************/
/* FIX OF taucs_available_memory_size()                             */
/*   returns size of memory available for allocation                */
/* (see taucs_memory.c)                                             */
/********************************************************************/

/* Redirect calls to avoid link error */
#define taucs_available_memory_size cgal_taucs_available_memory_size
#define taucs_malloc                malloc
#define taucs_free(ptr)             if (ptr != NULL) free(ptr)

/* LS 2007: if m_sys is meaningful, then we limit malloc test by 0.75*m_sys */
#define limit_memory(mem) ((mem) < m_max ? (mem) : m_max)

inline double taucs_available_memory_size() 
{
/* LS 2007: The original code below creates an infinite loop on Linux     */
/*          (optimistic memory allocation => malloc() never returns NULL) */
#ifndef OSTYPE_linux

  double m_sys;
  double m,m_low,m_high,m_tol;
  char*  p;
  double m_max;

  m_sys = taucs_system_memory_size();
  
  /* LS 2007: if m_sys is meaningful, then we limit malloc test by 0.75*m_sys */

  if (m_sys > 0) 
    m_max = floor(0.75 * m_sys); 
  else
    m_max = DBL_MAX;

  /* malloc test */

  m = 1048576.0;

  while ( (m < m_max-1) /* m_max not reached */
       && ((p=(char*) taucs_malloc( (size_t) limit_memory(m*2.0) )) != NULL) ) {
    taucs_printf("taucs_avail_memory_size: %.0lf Mb\n", limit_memory(m*2.0) / 1048576.0);
    taucs_free(p);
    m = limit_memory(m*2.0);
  }

  m_low  = m;
  m_high = limit_memory(m*2.0);
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

  return m;

#else /* LS 2007: workaround Linux optimistic memory allocation */

  double m_sys;   /* size of physical memory */
  double m;       /* size of memory available for allocation */

    m_sys  = (double) sysconf(_SC_PAGESIZE);
    m_sys *= (double) sysconf(_SC_PHYS_PAGES);

  /* we limit m by 0.75*m_sys */
  m = floor(0.75 * m_sys);

  taucs_printf((char*)"taucs_available_memory_size returns %lfMB\n", m/1048576.0);

  return m;
    
#endif
}


} // extern "C"


#endif // CGAL_TAUCS_FIX
