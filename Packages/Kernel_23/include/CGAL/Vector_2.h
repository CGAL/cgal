// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Andreas Fabri, Stefan Schirra

#ifndef CGAL_VECTOR_2_H
#define CGAL_VECTOR_2_H

#include <CGAL/Origin.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Vector_2 : public R_::Kernel_base::Vector_2
{
  typedef typename R_::RT             RT;
  typedef typename R_::Segment_2      Segment_2;
  typedef typename R_::Ray_2          Ray_2;
  typedef typename R_::Line_2         Line_2;
  typedef typename R_::Point_2        Point_2;
  typedef typename R_::Kernel_base::Vector_2  RVector_2;
public:
  typedef  R_                        R;

  Vector_2() {}

  Vector_2(const Point_2& a, const Point_2& b)
      : RVector_2(a, b) {}

  Vector_2(const RVector_2& v) : RVector_2(v) {}

  Vector_2(const Segment_2 &s) : RVector_2(s) {}

  Vector_2(const Ray_2 &r) : RVector_2(r) {}

  Vector_2(const Line_2 &l) : RVector_2(l) {}

  Vector_2(const Null_vector &v) : RVector_2(v) {}

  Vector_2(const RT &x, const RT &y) : RVector_2(x,y) {}

  Vector_2(const RT &x, const RT &y, const RT &w) : RVector_2(x,y,w) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_VECTOR_2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Vector_2<R> &v)
{
  typedef typename  R::Kernel_base::Vector_2  RVector_2;
  return os << static_cast<const RVector_2&>(v);
}
#endif // CGAL_NO_OSTREAM_INSERT_VECTOR_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_VECTOR_2
template < class R >
std::istream &
operator>>(std::istream &is, Vector_2<R> &p)
{
  typedef typename  R::Kernel_base::Vector_2  RVector_2;
  return is >> static_cast<RVector_2&>(p);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_VECTOR_2

CGAL_END_NAMESPACE

#endif // CGAL_VECTOR_2_H
