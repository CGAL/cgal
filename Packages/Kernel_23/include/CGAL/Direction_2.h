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
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Stefan Schirra
 
#ifndef CGAL_DIRECTION_2_H
#define CGAL_DIRECTION_2_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Direction_2 : public R_::Kernel_base::Direction_2
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Vector_2              Vector_2;
  typedef typename R_::Line_2                Line_2;
  typedef typename R_::Ray_2                 Ray_2;
  typedef typename R_::Segment_2             Segment_2;
  typedef typename R_::Kernel_base::Direction_2      RDirection_2;
public:
  typedef  R_   R;

  Direction_2() {}

  Direction_2(const RDirection_2& d)
    : RDirection_2(d) {}

  Direction_2(const Vector_2& v)
    : RDirection_2(v) {}

  Direction_2(const Line_2& l)
    : RDirection_2(l) {}

  Direction_2(const Ray_2& r)
    : RDirection_2(r) {}

  Direction_2(const Segment_2& s)
    : RDirection_2(s) {}

  Direction_2(const RT &x, const RT &y)
    :  RDirection_2(x,y) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_DIRECTION_2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Direction_2<R> &d)
{
  typedef typename  R::Kernel_base::Direction_2  RDirection_2;
  return os << static_cast<const RDirection_2&>(d);
}
#endif // CGAL_NO_OSTREAM_INSERT_DIRECTION_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_DIRECTION_2
template < class R >
std::istream &
operator>>(std::istream &is, Direction_2<R> &p)
{
  typedef typename  R::Kernel_base::Direction_2  RDirection_2;
  return is >> static_cast<RDirection_2&>(p);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_DIRECTION_2

CGAL_END_NAMESPACE

#endif // CGAL_DIRECTION_2_H
