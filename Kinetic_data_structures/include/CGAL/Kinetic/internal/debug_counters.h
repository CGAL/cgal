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

#include <CGAL/Kinetic/basic.h>

namespace CGAL { namespace Kinetic {
namespace internal {
  CGAL_EXPORT extern unsigned int zero_certificates__;
  CGAL_EXPORT extern unsigned int function_degeneracies__;
  CGAL_EXPORT extern unsigned int io_errors__;
  CGAL_EXPORT extern unsigned int audit_failures__;

  CGAL_EXPORT void write_debug_counters(std::ostream &out);
}
} } //namespace CGAL::Kinetic
