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
// Author(s)     : Stefan Schirra
 
#ifndef CGAL_ISO_CUBOID_3_H
#define CGAL_ISO_CUBOID_3_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Iso_cuboid_3 : public R_::Kernel_base::Iso_cuboid_3
{
  typedef typename R_::RT                 RT;
  typedef typename R_::Point_3            Point_3;
  typedef typename R_::Kernel_base::Iso_cuboid_3  RIso_cuboid_3;
public:
  typedef          R_                    R;

  Iso_cuboid_3() {}

  Iso_cuboid_3(const RIso_cuboid_3&  r)
      : RIso_cuboid_3(r) {}

  Iso_cuboid_3(const Point_3& p, const Point_3& q)
   : RIso_cuboid_3(p,q) {}

  Iso_cuboid_3(const Point_3 &left,   const Point_3 &right,
               const Point_3 &bottom, const Point_3 &top,
               const Point_3 &far_,   const Point_3 &close)
   : RIso_cuboid_3(left, right, bottom, top, far_, close) {}

  Iso_cuboid_3(const RT& min_hx, const RT& min_hy, const RT& min_hz,
               const RT& max_hx, const RT& max_hy, const RT& max_hz, 
               const RT& hw)
   : RIso_cuboid_3(min_hx, min_hy, min_hz, max_hx, max_hy, max_hz, hw) {}

  Iso_cuboid_3(const RT& min_hx, const RT& min_hy, const RT& min_hz,
               const RT& max_hx, const RT& max_hy, const RT& max_hz)
   : RIso_cuboid_3(min_hx, min_hy, min_hz, max_hx, max_hy, max_hz) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_ISO_CUBOID_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Iso_cuboid_3<R>& r)
{
  typedef typename  R::Kernel_base::Iso_cuboid_3  RIso_cuboid_3;
  return  os << (const RIso_cuboid_3& )r; }
#endif // CGAL_NO_OSTREAM_INSERT_ISO_CUBOID_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_ISO_CUBOID_3
template < class R >
std::istream&
operator>>(std::istream& is, Iso_cuboid_3<R>& r)
{
  typedef typename  R::Kernel_base::Iso_cuboid_3  RIso_cuboid_3;
  is >> (RIso_cuboid_3& )r;
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_ISO_CUBOID_3

CGAL_END_NAMESPACE

#endif // CGAL_ISO_CUBOID_3_H
