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
// Author(s)     : Andreas Fabri

#ifndef CGAL_RAY_2_H
#define CGAL_RAY_2_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Ray_2 : public R_::Kernel_base::Ray_2
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Point_2               Point_2;
  typedef typename R_::Direction_2           Direction_2;
  typedef typename R_::Vector_2              Vector_2;
  typedef typename R_::Line_2                Line_2;
  typedef typename R_::Kernel_base::Ray_2    RRay_2;
public:
  typedef  R_   R;

  Ray_2() {}

  Ray_2(const RRay_2& r)
    : RRay_2(r) {}

  Ray_2(const Point_2 &sp, const Point_2 &secondp)
    : RRay_2(sp, secondp) {}

  Ray_2(const Point_2 &sp, const Direction_2 &d)
    : RRay_2(sp, d) {}

  Ray_2(const Point_2 &sp, const Vector_2 &v)
    : RRay_2(sp, v) {}

  Ray_2(const Point_2 &sp, const Line_2 &l)
    : RRay_2(sp, l) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_RAY_2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Ray_2<R> &r)
{
  typedef typename  R::Kernel_base::Ray_2  RRay_2;
  return os << static_cast<const RRay_2&>(r);
}
#endif // CGAL_NO_OSTREAM_INSERT_RAY_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_RAY_2
template < class R >
std::istream &
operator>>(std::istream &is, Ray_2<R> &r)
{
  typedef typename  R::Kernel_base::Ray_2  RRay_2;
  return is >> static_cast<RRay_2&>(r);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_RAY_2

CGAL_END_NAMESPACE

#endif  // CGAL_RAY_2_H
