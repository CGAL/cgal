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

#ifndef CGAL_POLYNOMIAL_INTERNAL_EVEN_MULTIPLICITY_H
#define CGAL_POLYNOMIAL_INTERNAL_EVEN_MULTIPLICITY_H

#include <CGAL/Polynomial/basic.h>
#ifdef CGAL_USE_CORE
#include <CGAL/CORE_Expr.h>
#endif

namespace CGAL { namespace POLYNOMIAL { namespace internal {

//! Compute the sign after a root.
/*!
  This has specializations for Explicit_roots.
*/
template <class K>
class Is_even_multiplicity
{
public:
  Is_even_multiplicity(){  }

  typedef bool result_type;
  typedef typename K::Root argument_type;

  template <class T>
  result_type operator()(const T &v) const
  {
    return v.is_even_multiplicity();
  }

  bool operator()(double) const
  {
    return false;
  }
  bool operator()(Interval_nt) const
  {
    return false;
  }
#ifdef CGAL_USE_CORE
  bool operator()(::CORE::Expr) const
  {
    return false;
  }
#endif
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif
