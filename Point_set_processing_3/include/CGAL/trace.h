// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
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
// Author(s) : Laurent Saboret

#ifndef CGAL_TRACE_H
#define CGAL_TRACE_H

#include <stdio.h>
#include <iostream>
#include <fstream>


// Trace macros
// ------------

#ifdef DEBUG_TRACE
  #define CGAL_TRACE  printf
#else
  #define CGAL_TRACE  if (false) printf
#endif

#ifdef DEBUG_TRACE
  #define CGAL_TRACE_STREAM  std::cerr
#else
  #define CGAL_TRACE_STREAM  if (false) std::cerr
#endif


#endif // CGAL_TRACE_H
