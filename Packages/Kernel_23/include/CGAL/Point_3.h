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

#ifndef CGAL_POINT_3_H
#define CGAL_POINT_3_H

#include <CGAL/Origin.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Point_3 : public R_::Kernel_base::Point_3
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Vector_3              Vector_3;
  typedef typename R_::Kernel_base::Point_3  RPoint_3;
public:
  typedef          R_                       R;

  Point_3() {}

  Point_3(const Origin& o)
      : RPoint_3(o) {}

#if 1
  Point_3(const RPoint_3& p)
      : RPoint_3(p) {}
#endif

  Point_3(const RT& x, const RT& y, const RT& z)
    : RPoint_3(x, y, z) {}

  Point_3(const RT& hx, const RT& hy, const RT& hz, const RT& hw)
    : RPoint_3(hx, hy, hz, hw) {}
};

template <class R>
inline
bool
operator==(const Origin& o, const Point_3<R>& p)
{ return p == o; }

template <class R>
inline
bool
operator!=(const Origin& o, const Point_3<R>& p)
{ return p != o; }

#if 0 //ndef CGAL_NO_OSTREAM_INSERT_POINT_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Point_3<R>& p)
{
  typedef typename  R::Kernel_base::Point_3  RPoint_3;
  return os << static_cast<const RPoint_3&>(p);
}
#endif // CGAL_NO_OSTREAM_INSERT_POINT_3

#if 0 //ndef CGAL_NO_ISTREAM_EXTRACT_POINT_3
template < class R >
std::istream& operator>>(std::istream& is, Point_3<R>& p)
{
  typedef typename  R::Kernel_base::Point_3  RPoint_3;
  return is >> static_cast<RPoint_3&>(p);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_POINT_3

CGAL_END_NAMESPACE

#endif // CGAL_POINT_3_H
