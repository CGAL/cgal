// Copyright (c) 2005  Stanford University (USA).
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
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>


#include <CGAL/Tools/Log.h>
#include <CGAL/Kinetic/internal/debug_counters.h>
#include <iostream>
namespace CGAL {
Log::State Log::state_;
} //namespace CGAL
namespace CGAL { namespace Kinetic { namespace internal {

  unsigned int function_degeneracies__=0;
  unsigned int zero_certificates__=0;
  unsigned int io_errors__=0;
  unsigned int audit_failures__=0;

  void write_debug_counters(std::ostream &out) {
    out << "Degeneracies " << function_degeneracies__ << std::endl;
    out << "Zero functions " << zero_certificates__ << std::endl;
    if (io_errors__ != 0) out << "I/O errors " << io_errors__ << std::endl;
    if (audit_failures__ != 0) out << "Audit failures " << audit_failures__ << std::endl;
  }

} } } //namespace CGAL::Kinetic::internal
