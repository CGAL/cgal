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

#ifndef CGAL_TRIANGLE_3_H
#define CGAL_TRIANGLE_3_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Triangle_3 : public R_::Kernel_base::Triangle_3
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Point_3               Point_3;
  typedef typename R_::Kernel_base::Triangle_3  RTriangle_3;
public:
  typedef          R_                       R;

  Triangle_3() {}

  Triangle_3(const RTriangle_3& t)
      : RTriangle_3(t) {}

  Triangle_3(const Point_3& p, const Point_3& q, const Point_3& r)
    : RTriangle_3(p,q,r) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_TRIANGLE_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Triangle_3<R>& t)
{
  typedef typename  R::Kernel_base::Triangle_3  RTriangle_3;
  return os << static_cast<const RTriangle_3&>(t);
}
#endif // CGAL_NO_OSTREAM_INSERT_TRIANGLE_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_TRIANGLE_3
template < class R >
std::istream&
operator>>(std::istream& is, Triangle_3<R>& t)
{
  typedef typename  R::Kernel_base::Triangle_3  RTriangle_3;
  return is >> static_cast<RTriangle_3&>(t);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TRIANGLE_3

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGLE_3_H
