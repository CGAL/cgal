// Copyright (c) 2005  Stanford University (USA).
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
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_DEBUG_COUNTERS_H_
#define CGAL_DEBUG_COUNTERS_H_

#include <CGAL/Kinetic/basic.h>
#include <CGAL/atomic.h>

namespace CGAL { namespace Kinetic {
namespace internal {
  
#ifdef CGAL_NO_ATOMIC
  typedef unsigned int atomic_unsigned_int;
#else
  typedef CGAL::cpp11::atomic<unsigned int> atomic_unsigned_int;
#endif
  

  inline atomic_unsigned_int& get_static_function_degeneracies()
  { 
    static atomic_unsigned_int function_degeneracies__ ;
    return function_degeneracies__; 
  }
  inline atomic_unsigned_int& get_static_zero_certificates()
  {
    static atomic_unsigned_int zero_certificates__ ;
    return zero_certificates__; 
  }
  inline atomic_unsigned_int& get_static_io_errors()
  {
    static atomic_unsigned_int io_errors__ ;
    return io_errors__;
  }
  inline atomic_unsigned_int& get_static_audit_failures()
  {
    static atomic_unsigned_int audit_failures__ ;
    return audit_failures__;
  }




  CGAL_EXPORT void write_debug_counters(std::ostream &out);
}
} } //namespace CGAL::Kinetic

#ifdef CGAL_HEADER_ONLY
#include <CGAL/Kinetic/internal/debug_counters_impl.h>
#endif // CGAL_HEADER_ONLY

#endif // CGAL_DEBUG_COUNTERS_H_
