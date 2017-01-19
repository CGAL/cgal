// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s) : Laurent Saboret

#ifndef CGAL_TRACE_H
#define CGAL_TRACE_H

#include <CGAL/license/Point_set_processing_3.h>


#include <stdio.h>
#include <stdarg.h>
#include <iostream>
#include <fstream>

/// \cond SKIP_IN_MANUAL

// Trace utilities
// ---------------

// print_stderr() = printf-like function to print to stderr
inline void print_stderr(const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);
}

// CGAL_TRACE() = printf-like function to print to stderr
// if DEBUG_TRACE is defined (ignored otherwise)
#ifdef DEBUG_TRACE
  #define CGAL_TRACE  print_stderr
#else
  #define CGAL_TRACE  if (false) print_stderr
#endif

// CGAL_TRACE_STREAM = C++ stream that prints to std::cerr 
// if DEBUG_TRACE is defined (ignored otherwise)
#ifdef DEBUG_TRACE
  #define CGAL_TRACE_STREAM  std::cerr
#else
  #define CGAL_TRACE_STREAM  if (false) std::cerr
#endif

/// \endcond

#endif // CGAL_TRACE_H
