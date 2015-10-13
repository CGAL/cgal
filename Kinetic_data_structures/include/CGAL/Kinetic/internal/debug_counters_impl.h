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

#if defined(BOOST_MSVC)
#  pragma warning(disable:4251)
#endif
#include <iostream>

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

namespace CGAL { namespace Kinetic { namespace internal {

  CGAL_INLINE_FUNCTION
  void write_debug_counters(std::ostream &out) {
    out << "Degeneracies " << get_static_function_degeneracies() << std::endl;
    out << "Zero functions " << get_static_zero_certificates() << std::endl;
    if (get_static_io_errors() != 0) out << "I/O errors " << get_static_io_errors() << std::endl;
    if (get_static_audit_failures() != 0) out << "Audit failures " << get_static_audit_failures() << std::endl;
  }

} } } //namespace CGAL::Kinetic::internal
