// Copyright (c) 2000  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
//
// Author(s)     : Andreas Fabri, Herve Bronnimann

#ifndef CGAL_CARTESIAN_LINE_2_H
#define CGAL_CARTESIAN_LINE_2_H

#include <CGAL/array.h>
#include <CGAL/Handle_for.h>
#include <CGAL/predicates/kernel_ftC2.h>

namespace CGAL {

template < class R_ >
class LineC2
{
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Direction_2          Direction_2;
  typedef typename R_::Vector_2             Vector_2;
  typedef typename R_::Ray_2                Ray_2;
  typedef typename R_::Segment_2            Segment_2;
  typedef typename R_::Line_2               Line_2;

  typedef cpp11::array<FT, 3>               Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:

  typedef R_                                     R;

  LineC2() {}

  LineC2(const FT &a, const FT &b, const FT &c)
    : base(CGAL::make_array(a, b, c)) {}
  
  typename R_::Boolean   operator==(const LineC2 &l) const;
  typename R_::Boolean   operator!=(const LineC2 &l) const;

  const FT & a() const
  {
      return get_pointee_or_identity(base)[0];
  }
  const FT & b() const
  {
      return get_pointee_or_identity(base)[1];
  }
  const FT & c() const
  {
      return get_pointee_or_identity(base)[2];
  }

};

template < class R >
CGAL_KERNEL_INLINE
typename R::Boolean
LineC2<R>::operator==(const LineC2<R> &l) const
{
  if (CGAL::identical(base, l.base))
      return true;
  return equal_lineC2(a(), b(), c(), l.a(), l.b(), l.c());
}

template < class R >
inline
typename R::Boolean
LineC2<R>::operator!=(const LineC2<R> &l) const
{
  return ! (*this == l);
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_LINE_2_H
