// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
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
//                 Sven Schoenherr

#ifndef CGAL_CIRCLE_2_H
#define CGAL_CIRCLE_2_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Circle_2 : public R_::Kernel_base::Circle_2
{
  typedef typename R_::FT                    FT;
  typedef typename R_::Point_2               Point_2;
  typedef typename R_::Kernel_base::Circle_2  RCircle_2;
public:
  typedef  R_   R;

    Circle_2() {}

    Circle_2(const RCircle_2& t)
      : RCircle_2(t) {}

    Circle_2(const Point_2 &center, const FT &squared_radius,
             const Orientation &orientation)
      : RCircle_2(center, squared_radius, orientation) {}

    Circle_2(const Point_2 &center, const FT &squared_radius)
      : RCircle_2(center, squared_radius, COUNTERCLOCKWISE) {}

    Circle_2(const Point_2 &p, const Point_2 &q, const Point_2 &r)
      : RCircle_2(p,q,r) {}

    Circle_2(const Point_2 & p, const Point_2 & q,
             const Orientation &orientation)
      : RCircle_2(p,q,orientation) {}

    Circle_2(const Point_2 & p, const Point_2 & q)
      : RCircle_2(p,q,COUNTERCLOCKWISE) {}

    Circle_2(const Point_2 & center, const Orientation& orientation)
      : RCircle_2(center,FT(0),orientation) {}

    Circle_2(const Point_2 & center)
      : RCircle_2(center,FT(0),COUNTERCLOCKWISE) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_CIRCLE_2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Circle_2<R> &c)
{
  typedef typename R::Kernel_base::Circle_2  RCircle_2;
  return os << (const RCircle_2&)c;
}

#endif // CGAL_NO_OSTREAM_INSERT_CIRCLE_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_CIRCLE_2
template < class R >
std::istream &
operator>>(std::istream &is, Circle_2<R> &c)
{
  typedef typename R::Kernel_base::Circle_2  RCircle_2;
  return is >> (RCircle_2&)c;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_CIRCLE_2

CGAL_END_NAMESPACE

#endif  // CGAL_CIRCLE_2_H
