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
/*                                                                  */
/* FIX OF TAUCS BUGGY FUNCTION taucs_available_memory_size()        */
/*                                                                  */
/********************************************************************/

/* taucs_available_memory_size/cgal_taucs_available_memory_size     */
/*   returns size of memory available for allocation                */
inline double cgal_taucs_available_memory_size()
{
    double m;       /* size of memory available for allocation */

#ifdef OSTYPE_linux
    /* taucs_available_memory_size() is buggy on Linux 2.6 */
    /* It returns only 1% of the actual memory             */
    double m_sys;   /* size of physical memory */
    m_sys  = (double) sysconf(_SC_PAGESIZE);
    m_sys *= (double) sysconf(_SC_PHYS_PAGES);

    /* we limit m by 0.75*m_sys */
    m = floor(0.75 * m_sys);
#else
    m = taucs_available_memory_size();
#endif

    taucs_printf("cgal_taucs_available_memory_size returns %lfMB\n", m/1048576.0);

    return m;
}

/* Redirect calls to cgal_taucs_available_memory_size()             */
/* to taucs_available_memory_size()                                 */
#define taucs_available_memory_size cgal_taucs_available_memory_size

} // extern "C"


#endif // CGAL_TAUCS_FIX
