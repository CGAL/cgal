// Copyright (c) 2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Laurent Saboret
//                 (based on Memory_sizer.h by Sylvain Pion and Andreas Fabri)

#ifndef CGAL_PEAK_MEMORY_SIZER_H
#define CGAL_PEAK_MEMORY_SIZER_H

#include <CGAL/trace.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/Taucs_fix.h>
#include <CGAL/surface_reconstruction_points_assertions.h>

#include <deque>
#include <cmath>
#include <cfloat>
#include <climits>

namespace CGAL {


/// Peak_memory_sizer extends Memory_sizer with new memory statistics.
struct Peak_memory_sizer : public Memory_sizer
{
    typedef std::size_t   size_type;

    /// Gets the peak memory used by the process.
    /// Both the virtual memory size and the resident size.
    size_type peak_virtual_size()  const { return get_peak_memory(true); }
    size_type peak_resident_size() const { return get_peak_memory(false); }

    /// Gets the number of large free memory blocks.
    size_type count_free_memory_blocks(size_type block_size) const
    {
      double m_sys; /* physical system memory size */
      double m; /* allocated memory */
      double m_max; /* malloc test limit */

      /* get physical system memory size */
      m_sys = taucs_system_memory_size();

      /* LS 2007: if m_sys is meaningful, then we limit malloc test by 0.75*m_sys */
      /*          to avoid an infinite loop on Linux (malloc() never returns NULL */
      /*          due to "optimistic memory allocation")                          */
      if (m_sys > 0)
        m_max = floor(0.75 * m_sys);
      else
        m_max = DBL_MAX;

      // Allocate all memory blocks >= block_size
      // while keeping the total memory allocated <= m_max.
      m = 0;
      std::deque<void*> blocks;
      void* block;
      while ( (m + block_size <= m_max) /* m_max not reached */
        && ((block = malloc(block_size)) != NULL) )
      {
        m += block_size;
        //CGAL_TRACE("allocated large memory blocks up to %.0lf Mb\n", m / 1048576.0);
        blocks.push_back(block);
      }

      // Return value
      size_type count = blocks.size();

        // Free large memory blocks
        for (size_type i=0; i<count; i++)
          free(blocks[i]);

      //CGAL_TRACE("%ld large blocks are free\n", count);
      return count;
    }

    /// Give size of largest block available for allocation.
    // (based on taucs_available_memory_size() by S. Toledo)
    size_t largest_free_block() const
    {
#ifdef _WIN64
  // TEMPORARY HACK to avoid an infinite loop in malloc(2 Gb) in Poisson MFC demo on Windows XP 64 bits!
  // Note: the code works fine under the debugger!
  return 0;
#endif

      double m_sys; /* physical system memory size */
      double m, /* allocated memory */
             m_low,m_high,m_tol; /* memory range */
      char*  p;
      double m_max; /* malloc test limit */

      /* get physical system memory size */
      m_sys = taucs_system_memory_size();

      /* LS 2007: if m_sys is meaningful, then we limit malloc test by 0.75*m_sys */
      /*          to avoid an infinite loop on Linux (malloc() never returns NULL */
      /*          due to "optimistic memory allocation")                          */
      if (m_sys > 0)
        m_max = floor(0.75 * m_sys);
      else
        m_max = DBL_MAX;

      /* malloc test */

      m = 1048576.0;

      while ( (m < m_max) /* m_max not reached */
           && ((p=(char*) malloc( (size_t) (std::min)(m_max,m*2.0) )) != NULL) ) {
         //CGAL_TRACE("largest_free_block: %.0lf Mb\n", (std::min)(m_max,m*2.0) / 1048576.0);
        free(p);
        m = (std::min)(m_max,m*2.0);
      }

      m_low  = m;
      m_high = (std::min)(m_max,m*2.0);
      m_tol  = m / 128.0;

      while ( m_high - m_low > m_tol ) {
        m = m_low + ( (m_high-m_low)/2.0 );
        //CGAL_TRACE("largest_free_block: [%.0lf %.0lf %.0lf]\n",
        //  	       m_low  / 1048576.0,
        //  	       m      / 1048576.0,
        //  	       m_high / 1048576.0);
        if ( (p=(char*) malloc( (size_t) m )) != NULL )
          m_low = m;
        else
          m_high = m;
        free(p);
      }

      m = m_low;

      //CGAL_TRACE("largest_free_block: malloc test=%.0lf MB max test=%.0lf MB\n",
      //    	     m / 1048576.0,
      //    	     m_max / 1048576.0);
      return (size_t) m;
    }

private:

  size_type get_peak_memory (bool virtual_size)  const
  {
#ifdef _MSC_VER
    DWORD pid = GetCurrentProcessId();
    size_type result;
    HANDLE hProcess;
    PROCESS_MEMORY_COUNTERS pmc;
    hProcess = OpenProcess(  PROCESS_QUERY_INFORMATION |
                             PROCESS_VM_READ,
                             FALSE, pid );
    if ( GetProcessMemoryInfo( hProcess, &pmc, sizeof(pmc)) )
    {
//CGAL_TRACE("    Peak_memory_sizer: WorkingSetSize=%ld Mb\n",              pmc.WorkingSetSize>>20);
//CGAL_TRACE("    Peak_memory_sizer: PagefileUsage=%ld Mb\n",               pmc.PagefileUsage>>20);
//CGAL_TRACE("    Peak_memory_sizer: PeakWorkingSetSize=%ld Mb\n",          pmc.PeakWorkingSetSize>>20);
//CGAL_TRACE("    Peak_memory_sizer: PeakPagefileUsage=%ld Mb\n",           pmc.PeakPagefileUsage>>20);

      // LS 10/2008: PeakPagefileUsage seems unreliable, thus we use an approximation:
      size_t memory_paged_out = (pmc.PagefileUsage>pmc.WorkingSetSize) ? (pmc.PagefileUsage-pmc.WorkingSetSize) : 0;
      size_t approximate_peak_virtual_size = pmc.PeakWorkingSetSize + memory_paged_out;
//CGAL_TRACE("    Peak_memory_sizer: approximate_peak_virtual_size=%ld Mb\n", approximate_peak_virtual_size>>20);

      result = virtual_size ? approximate_peak_virtual_size : pmc.PeakWorkingSetSize;
    }

    CloseHandle( hProcess );
    return result;

#else
    // Not yet implemented
    return 0;
#endif
  }
};


} //namespace CGAL

#endif // CGAL_PEAK_MEMORY_SIZER_H
