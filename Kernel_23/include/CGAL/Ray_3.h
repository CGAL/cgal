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
 
#ifndef CGAL_RAY_3_H
#define CGAL_RAY_3_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Ray_3 : public R_::Kernel_base::Ray_3
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Point_3               Point_3;
  typedef typename R_::Direction_3           Direction_3;
  typedef typename R_::Vector_3              Vector_3;
  typedef typename R_::Line_3                Line_3;
  typedef typename R_::Kernel_base::Ray_3    RRay_3;
public:
  typedef          R_                       R;

  Ray_3() {}

  Ray_3(const RRay_3& r)
      : RRay_3(r) {}

  Ray_3(const Point_3& sp, const Point_3& secondp)
    : RRay_3(sp, secondp) {}

  Ray_3(const Point_3& sp, const Vector_3& v)
    : RRay_3(sp, v) {}

  Ray_3(const Point_3& sp, const Direction_3& d)
    : RRay_3(sp, d) {}

  Ray_3(const Point_3& sp, const Line_3& l)
    : RRay_3(sp, l) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_RAY_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Ray_3<R>& r)
{
  typedef typename  R::Kernel_base::Ray_3  RRay_3;
  return os << static_cast<const RRay_3&>(r);
}
#endif // CGAL_NO_OSTREAM_INSERT_RAY_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_RAY_3
template < class R >
std::istream&
operator>>(std::istream& is, Ray_3<R>& r)
{
  typedef typename  R::Kernel_base::Ray_3  RRay_3;
  return is >> static_cast<RRay_3&>(r);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_RAY_3

CGAL_END_NAMESPACE

#endif // CGAL_RAY_3_H
