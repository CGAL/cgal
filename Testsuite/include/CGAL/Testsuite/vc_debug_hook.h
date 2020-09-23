// Copyright (c) 2008 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola
//
// This is used by the testsuite to prevent Visual C++ from poping up an error window.
//

#ifndef CGAL_VC_DEBUG_HOOK_H
#define CGAL_VC_DEBUG_HOOK_H

#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <crtdbg.h>


namespace
{
  int CGAL_report_hook ( int type, char* msg, int* retval )
  {
    if ( type == _CRT_ASSERT )
    {
      std::fprintf(stderr,msg);
      *retval = 0 ;
      std::exit(255);
    }
    return 1 ;
  }

  void CGAL_handle_signal( int n )
  {
    switch(n)
    {
      case SIGSEGV: std::fprintf(stderr,"In CGAL_handle_signal, Program received signal SIGSEGV: Segmentation Fault."); break ;
      case SIGFPE : std::fprintf(stderr,"In CGAL_handle_signal, Program received signal SIGFPE: Floating Point Execption."); break ;
      case SIGILL : std::fprintf(stderr,"In CGAL_handle_signal, Program received signal SIGILL: Illegal Instruction."); break ;
      default:
        std::fprintf(stderr,"In CGAL_handle_signal, Program received signal %d", n); break ;
    }

    std::exit(128+n);
  }

  struct CGAL_DebugHook
  {
    CGAL_DebugHook()
    {
      _CrtSetReportHook(CGAL_report_hook);

      // This is OK for unattended runs but will prevent the IDE for trapping the signal
      std::signal(SIGSEGV,CGAL_handle_signal);
      std::signal(SIGFPE ,CGAL_handle_signal);
      std::signal(SIGILL ,CGAL_handle_signal);
    }
  } ___ ;
}

#endif // CGAL_VC_DEBUG_HOOK_H
