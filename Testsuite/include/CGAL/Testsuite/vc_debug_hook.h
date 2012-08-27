// Copyright (c) 2008 GeometryFactory (France).
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
      case SIGSEGV: std::fprintf(stderr,"Program recieved signal SIGSEGV: Segmentation Fault."); break ;
      case SIGFPE : std::fprintf(stderr,"Program recieved signal SIGFPE: Floating Point Execption."); break ;
      case SIGILL : std::fprintf(stderr,"Program recieved signal SIGILL: Illegal Instruction."); break ;
      default:
        std::fprintf(stderr,"Program recieved signal %d", n); break ;
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
