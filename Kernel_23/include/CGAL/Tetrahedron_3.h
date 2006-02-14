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

#ifndef CGAL_TETRAHEDRON_3_H
#define CGAL_TETRAHEDRON_3_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Tetrahedron_3 : public R_::Kernel_base::Tetrahedron_3
{
  typedef typename R_::Point_3             Point_3;
  typedef typename R_::Kernel_base::Tetrahedron_3  RTetrahedron_3;
public:
  typedef          R_                       R;

  Tetrahedron_3() {}

  Tetrahedron_3(const RTetrahedron_3& t)
      : RTetrahedron_3(t) {}

  Tetrahedron_3(const Point_3& p,
                const Point_3& q,
                const Point_3& r,
                const Point_3& s)
    : RTetrahedron_3(p,q,r,s) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_TETRAHEDRON_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Tetrahedron_3<R>& t)
{
  typedef typename  R::Kernel_base::Tetrahedron_3  RTetrahedron_3;
  return os << static_cast<const RTetrahedron_3&>(t);
}
#endif // CGAL_NO_OSTREAM_INSERT_TETRAHEDRON_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_TETRAHEDRON_3
template < class R >
std::istream&
operator>>(std::istream& is, Tetrahedron_3<R>& t)
{
  typedef typename  R::Kernel_base::Tetrahedron_3  RTetrahedron_3;
  return is >> static_cast<RTetrahedron_3&>(t);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TETRAHEDRON_3

CGAL_END_NAMESPACE

#endif  // CGAL_TETRAHEDRON_3_H
