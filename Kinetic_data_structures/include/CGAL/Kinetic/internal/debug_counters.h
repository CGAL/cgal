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

namespace CGAL { namespace Kinetic {
namespace internal {
  
#ifdef CGAL_HEADER_ONLY
  
  inline unsigned int& get_static_function_degeneracies()
  { 
    static unsigned int function_degeneracies__ = 0;
    return function_degeneracies__; 
  }
  inline unsigned int& get_static_zero_certificates()
  {
    static unsigned int zero_certificates__ = 0;
    return zero_certificates__; 
  }
  inline unsigned int& get_static_io_errors()
  {
    static unsigned int io_errors__ = 0;
    return io_errors__;
  }
  inline unsigned int& get_static_audit_failures()
  {
    static unsigned int audit_failures__ = 0;
    return audit_failures__;
  }

#else // CGAL_HEADER_ONLY
  
  CGAL_EXPORT extern unsigned int function_degeneracies__;
  CGAL_EXPORT extern unsigned int zero_certificates__;
  CGAL_EXPORT extern unsigned int io_errors__;
  CGAL_EXPORT extern unsigned int audit_failures__;

  inline unsigned int& get_static_function_degeneracies()
  { return function_degeneracies__; }
  inline unsigned int& get_static_zero_certificates()
  { return zero_certificates__; }
  inline unsigned int& get_static_io_errors()
  { return io_errors__; }
  inline unsigned int& get_static_audit_failures()
  { return audit_failures__; }

#endif // CGAL_HEADER_ONLY


  CGAL_EXPORT void write_debug_counters(std::ostream &out);
}
} } //namespace CGAL::Kinetic

#ifdef CGAL_HEADER_ONLY
#include <CGAL/Kinetic/internal/debug_counters_impl.h>
#endif // CGAL_HEADER_ONLY

#endif // CGAL_DEBUG_COUNTERS_H_
