// Copyright (c) 2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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

// This has only been implemented for MacOSX/Darwin, Linux and VC++ for now.
#if !defined _MSC_VER && !defined __linux__ && !defined __APPLE__

#include <iostream>

namespace CGAL {

struct Memory_sizer
{
    typedef std::size_t   size_type;
    size_type virtual_size()  const { return 0; }
    size_type resident_size() const { return 0; }
};

} //namespace CGAL

#else // defined _MSC_VER ||  defined __linux__ || defined __APPLE__

#if defined _MSC_VER
#  include <windows.h>
#  include "psapi.h"

// auto-link with psapi.lib
#  define CGAL_LIB_NAME psapi
#  define CGAL_AUTO_LINK_NOMANGLE
#  include <CGAL/auto_link/auto_link.h>

#elif defined __linux__ 
#  include <fstream>
#  include <cstddef>
#  include <unistd.h>
#elif defined __APPLE__
#include <mach/task.h>
#include <mach/mach_init.h>
#endif

namespace CGAL {

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
      // PagefileUsage seems not very precise, thus check it against WorkingSetSize:
      size_t approximate_virtual_size = (std::max)(pmc.PagefileUsage, pmc.WorkingSetSize);
      
      result = (virtual_size)? approximate_virtual_size : pmc.WorkingSetSize;
    }

    CloseHandle( hProcess );
    return result;

#elif defined __linux__
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
    CGAL_assertion(!f.bad());

    f >> pid >> name >> state >> ppid >> pgrp >> session >> tty >> tpgid >> flags;
    f >> minflt >> cminflt >> majflt >> cmajflt >> utime >> stime >> cutime;
    f >> cstime >> counter >> priority >> timeout >> itrealvalue >> starttime;
    f >> vsize >> rss;

    return virtual_size ? vsize : rss * getpagesize();

#else // __APPLE__ is defined

    // http://miknight.blogspot.com/2005/11/resident-set-size-in-mac-os-x.html
		// This is highly experimental. But still better than returning 0.
		// It appears that we might need certain 'rights' to get access to the kernel
		// task... It works if you have admin rights apparently
		// (though non-root of course!). I haven't tested with non-admin user.
    // -- Samuel Hornus

    task_t task = MACH_PORT_NULL;
		// The task_for_pid() seems to be time consuming (looking at the source
		// in xnu-source/bsd/vm/vm_unix.c
		// TODO: so it may be a good idea to cache the resulting 'task'
    if (task_for_pid(current_task(), getpid(), &task) != KERN_SUCCESS)
        return 0;
		// It seems to me that after calling :
		// task_for_pid(current_task(), getpid(), &task)
		// we should have (task == current_task())
		// ==> TODO: Check if this is indeed the case
      
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
      
    task_info(task, TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count);
		// TODO: test if the following line works...
    //task_info(current_task(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count);
#if 0
    std::cerr << "PAGE SIZE IS " << getpagesize() << std::endl
    << " RESIDENT SIZE IS " << t_info.resident_size << std::endl
    << " VIRTUAL SIZE IS " << t_info.virtual_size << std::endl;
#endif
    return virtual_size ? t_info.virtual_size : t_info.resident_size;
#endif
  }
};

} //namespace CGAL

#endif

#endif
