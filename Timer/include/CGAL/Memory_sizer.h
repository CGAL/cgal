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
// Author(s)     : Sylvain Pion, Andreas Fabri

#ifndef CGAL_MEMORY_SIZER_H
#define CGAL_MEMORY_SIZER_H

#include <CGAL/basic.h>

// This has only been implemented for Linux and VC++ for now.
#if !defined _MSC_VER && !defined __linux__

#include <iostream>

CGAL_BEGIN_NAMESPACE

struct Memory_sizer
{
    typedef std::size_t   size_type;
    size_type virtual_size()  const { return 0; }
    size_type resident_size() const { return 0; }
};

CGAL_END_NAMESPACE

#else // defined _MSC_VER ||  defined __linux__

#if defined _MSC_VER
#  include <windows.h>
#  include "psapi.h"
#elif defined __linux__ 
#  include <fstream>
#  include <cstddef>
#  include <unistd.h>
#endif

CGAL_BEGIN_NAMESPACE

// A class giving access to the memory currently used by the process.
// Both the virtual memory size and the resident size.
// I put it in a class instead of free functions for similarity with Timer,
// and in case we want to store some state.

struct Memory_sizer
{
    typedef std::size_t   size_type;

    size_type virtual_size()  const { return get(true); }
    size_type resident_size() const { return get(false); }

private:

  size_type get (bool virtual_size)  const 
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
      result = (virtual_size)? pmc.PagefileUsage : pmc.WorkingSetSize;
    }

    CloseHandle( hProcess );
    return result;

#else
    // Extract of "man proc" under Linux :
    //
    //            vsize %u Virtual memory size
    //
    //            rss %u Resident Set Size: number of pages
    //                   the process has in real memory,
    //                   minus 3 for administrative purposes.
    //                   This is just the pages which count
    //                   towards text, data, or stack space.
    //                   This does not include pages which
    //                   have not been demand-loaded in, or
    //                   which are swapped out.
    //
    // Note : the following may be buggy in case of space in the executable name...

    int pid;
    char name[1024];
    char state;
    int ppid, pgrp, session, tty, tpgid;
    unsigned flags, minflt, cminflt, majflt, cmajflt;
    int utime, stime, cutime, cstime, counter, priority, timeout;
    unsigned itrealvalue, starttime;
    size_type vsize = 0, rss = 0;

    std::ifstream f("/proc/self/stat");
    CGAL_assertion(f);

    f >> pid >> name >> state >> ppid >> pgrp >> session >> tty >> tpgid >> flags;
    f >> minflt >> cminflt >> majflt >> cmajflt >> utime >> stime >> cutime;
    f >> cstime >> counter >> priority >> timeout >> itrealvalue >> starttime;
    f >> vsize >> rss;

    return virtual_size ? vsize : rss * getpagesize();
#endif
  }
};

CGAL_END_NAMESPACE

#endif

#endif
