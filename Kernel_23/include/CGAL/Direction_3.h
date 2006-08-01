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
// $URL$
// $Id$
// 
//
// Author(s)     : Andreas Fabri, Stefan Schirra
 
#ifndef CGAL_DIRECTION_3_H
#define CGAL_DIRECTION_3_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Direction_3 : public R_::Kernel_base::Direction_3
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Vector_3              Vector_3;
  typedef typename R_::Line_3                Line_3;
  typedef typename R_::Ray_3                 Ray_3;
  typedef typename R_::Segment_3             Segment_3;
  typedef typename R_::Aff_transformation_3  Aff_transformation_3;
  typedef typename R_::Kernel_base::Direction_3 RDirection_3;
public:

  typedef RDirection_3 Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef          R_                       R;

  Direction_3() {}

  Direction_3(const RDirection_3& d)
    : RDirection_3(d) {}

  Direction_3(const Vector_3& v)
    : RDirection_3(v) {}

  Direction_3(const Line_3& l)
    : RDirection_3(l) {}

  Direction_3(const Ray_3& r)
    : RDirection_3(r) {}

  Direction_3(const Segment_3& s)
    : RDirection_3(s) {}

  Direction_3(const RT& hx, const RT& hy, const RT& hz)
    : RDirection_3(hx, hy, hz) {}

  Direction_3 transform(const Aff_transformation_3 &t) const
  {
    return t.transform(*this);
  }
 
  Direction_3
  operator-() const
  {
    return R().construct_opposite_direction_3_object()(*this);
  } 
  

};

#ifndef CGAL_NO_OSTREAM_INSERT_DIRECTION_3
template < class R >
std::ostream& operator<<(std::ostream& os, const Direction_3<R>& d)
{
  typedef typename  R::Kernel_base::Direction_3 RDirection_3;
  return os << static_cast<const RDirection_3&>(d); }
#endif // CGAL_NO_OSTREAM_INSERT_DIRECTION_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_DIRECTION_3
template < class R >
std::istream& operator>>(std::istream& is, Direction_3<R>& p)
{
  typedef typename  R::Kernel_base::Direction_3  RDirection_3;
  return is >> static_cast<RDirection_3&>(p); }
#endif // CGAL_NO_ISTREAM_EXTRACT_DIRECTION_3

CGAL_END_NAMESPACE

#endif // CGAL_DIRECTION_3_H
