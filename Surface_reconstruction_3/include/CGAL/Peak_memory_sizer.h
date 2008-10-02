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

#include <CGAL/Memory_sizer.h>

CGAL_BEGIN_NAMESPACE

// Peak_memory_sizer extends Memory_sizer by giving access to the peak memory used by the process.
// Both the virtual memory size and the resident size.
struct Peak_memory_sizer : public Memory_sizer
{
    typedef std::size_t   size_type;

    size_type peak_virtual_size()  const { return get_peak_memory(true); }
    size_type peak_resident_size() const { return get_peak_memory(false); }

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

CGAL_END_NAMESPACE

#endif // CGAL_PEAK_MEMORY_SIZER_H
