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
// Author(s)     : Andreas Fabri, Stefan Schirra
 
#ifndef CGAL_PLANE_3_H
#define CGAL_PLANE_3_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Plane_3 : public R_::Kernel_base::Plane_3
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Point_3               Point_3;
  typedef typename R_::Direction_3           Direction_3;
  typedef typename R_::Vector_3              Vector_3;
  typedef typename R_::Segment_3             Segment_3;
  typedef typename R_::Line_3                Line_3;
  typedef typename R_::Ray_3                 Ray_3;
  typedef typename R_::Kernel_base::Plane_3  RPlane_3;
public:
  typedef          R_                       R;

  Plane_3() {}

  Plane_3(const RPlane_3& p)
      : RPlane_3(p) {}

  Plane_3(const Point_3& p, const Point_3& q, const Point_3& r)
    : RPlane_3(p,q,r) {}

  Plane_3(const Point_3& p, const Direction_3& d)
    : RPlane_3(p,d) {}

  Plane_3(const Point_3& p, const Vector_3& v)
    : RPlane_3(p,v) {}

  Plane_3(const RT& a, const RT& b, const RT& c, const RT& d)
    : RPlane_3(a,b,c,d) {}

  Plane_3(const Line_3& l, const Point_3& p)
    : RPlane_3(l,p) {}

  Plane_3(const Segment_3& s, const Point_3& p)
    : RPlane_3(s,p) {}

  Plane_3(const Ray_3& r, const Point_3& p)
    : RPlane_3(r,p) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_PLANE_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Plane_3<R>& p)
{
  typedef typename  R::Kernel_base::Plane_3  RPlane_3;
  return os << static_cast<const RPlane_3&>(p);
}
#endif // CGAL_NO_OSTREAM_INSERT_PLANE_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_PLANE_3
template < class R >
std::istream&
operator>>(std::istream& is, Plane_3<R>& t)
{
  typedef typename  R::Kernel_base::Plane_3  RPlane_3;
  return is >> static_cast<RPlane_3&>(t);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_PLANE_3

CGAL_END_NAMESPACE

#endif // CGAL_PLANE_3_H
