// Copyright (c) 2005-2008  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Sylvain Pion, Mariette Yvinec, Laurent Rineau

#ifndef CGAL_ROBUST_CONSTRUCTION_H
#define CGAL_ROBUST_CONSTRUCTION_H

#include <CGAL/config.h>

namespace CGAL {

// This template class is a functor adaptor targetting geometric constructions.
//
// They are "robust" in the following sense : the input and output are
// approximate (doubles), but the internal computation tries to guarantees the
// approximation quantitatively.
//
// This is especially useful in order to minimize the bad things that happen
// in close-to-degenerate cases (e.g. computing the circumcenter of an almost
// flat tetrahedron).
//
// The implementation strategy is to:
// - convert the input to Lazy kernel objects
// - perform the computation with it
// - convert back to double in a guaranteed way (since Lazy_exact_nt guarantees
//   a precision on to_double).

// TODO :
// - possible improvement by avoiding the constructions of Lazy_exact_nt's.

template <class EC, class A2E, class E2A, class Result_type>
class Robust_construction
{
  EC ec;
  A2E a2e;
  E2A e2a;

public:

  typedef Result_type  result_type;

  typedef EC    Exact_construction;
  typedef A2E   Approximate_to_exact_converter;
  typedef E2A   Exact_to_approximate_converter;

  template <typename... Args>
  result_type
  operator()(const Args&... args) const
  { return e2a(ec(a2e(args)...)); }
};

} //namespace CGAL

#endif // CGAL_ROBUST_CONSTRUCTION_H
