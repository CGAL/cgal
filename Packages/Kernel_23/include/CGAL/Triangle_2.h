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
// Author(s)     : Andreas Fabri
 

#ifndef CGAL_TRIANGLE_2_H
#define CGAL_TRIANGLE_2_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Triangle_2 : public R_::Kernel_base::Triangle_2
{
  typedef typename R_::Point_2          Point_2;
  typedef typename R_::Kernel_base::Triangle_2  RTriangle_2;
public:
  typedef  R_                          R;

  Triangle_2() {}

  Triangle_2(const RTriangle_2& t)
      : RTriangle_2(t) {}

  Triangle_2(const Point_2 &p, const Point_2 &q, const Point_2 &r)
      : RTriangle_2(p,q,r) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_TRIANGLE_2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Triangle_2<R> &t)
{
  typedef typename  R::Kernel_base::Triangle_2  RTriangle_2;
  return os << (const RTriangle_2&)t;
}
#endif // CGAL_NO_OSTREAM_INSERT_TRIANGLE_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_TRIANGLE_2
template < class R >
std::istream &
operator>>(std::istream &is, Triangle_2<R> &t)
{
  typedef typename  R::Kernel_base::Triangle_2  RTriangle_2;
  return is >> (RTriangle_2&)t;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TRIANGLE_2

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGLE_2_H
