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

#ifndef CGAL_VECTOR_3_H
#define CGAL_VECTOR_3_H

#include <CGAL/Origin.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Vector_3 : public R_::Kernel_base::Vector_3
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Segment_3             Segment_3;
  typedef typename R_::Ray_3                 Ray_3;
  typedef typename R_::Line_3                Line_3;
  typedef typename R_::Point_3               Point_3;
  typedef typename R_::Kernel_base::Vector_3         RVector_3;
public:
  typedef          R_                       R;

  Vector_3() {}

  Vector_3(const Point_3& a, const Point_3& b)
    : RVector_3(a, b) {}

  Vector_3(const Segment_3& s)
      : RVector_3(s) {}

  Vector_3(const Ray_3& r)
      : RVector_3(r) {}

  Vector_3(const Line_3& l)
      : RVector_3(l) {}

  Vector_3(const RVector_3& v)
      : RVector_3(v) {}

  Vector_3(const Null_vector& v)
      : RVector_3(v) {}

  Vector_3(const RT& x, const RT& y, const RT& z)
    : RVector_3(x, y, z) {}

  Vector_3(const RT& x, const RT& y, const RT& z, const RT& w)
    : RVector_3(x, y, z, w) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_VECTOR_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Vector_3<R>& v)
{
  typedef typename  R::Kernel_base::Vector_3  RVector_3;
  return os << static_cast<const RVector_3&>(v);
}
#endif // CGAL_NO_OSTREAM_INSERT_VECTOR_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_VECTOR_3
template < class R >
std::istream&
operator>>(std::istream& is, Vector_3<R>& p)
{
  typedef typename  R::Kernel_base::Vector_3  RVector_3;
  return is >> static_cast<RVector_3&>(p);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_VECTOR_3

CGAL_END_NAMESPACE

#endif // CGAL_VECTOR_3_H
